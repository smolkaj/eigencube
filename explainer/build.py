"""Renders the chapters in parallel, then stitches them into one captioned film.

    python build.py            # 1080p30 final cut -> build/eigencube-explainer.mp4
    python build.py --draft    # fast 480p15 preview -> build/draft.mp4

Only chapters whose inputs changed since their last render are rendered again.
"""

import argparse
import ast
import hashlib
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"
sys.path.insert(0, str(HERE))
from kit import duration, speakable  # noqa: E402
from scenes import SCENES  # noqa: E402

CHAPTERS = [scene.__name__ for scene in SCENES]
CAPTION_WIDTH = 50  # Characters, so that a caption fits on one line in common players.
FINAL_FPS = 30
MIN_CAPTION_SECONDS = 0.8
MIN_CAPTION = 24  # Characters; shorter sentences are joined to a neighbor, so they don't flash by.


def render(scene, flags, media):
    log = media / f"{scene}.log"
    # Caching stays off: after replaying a cached animation, Manim silently drops any sound added
    # before the next one, which cost three narration lines their audio.
    with log.open("w") as out:
        rendered = subprocess.run([sys.executable, "-m", "manim", "render", *flags,
                                   "--disable_caching", "--media_dir", str(media),
                                   str(HERE / "scenes.py"), scene], stdout=out, stderr=out)
    if rendered.returncode:
        sys.exit(f"{scene} failed to render; see {log}")


def encode(video, encoded):
    """A chapter's picture as it goes into the film. Encoded chapter by chapter, as each finishes
    rendering and in parallel with the others, the film's picture is then a plain join."""
    ffmpeg_to(encoded, ["-i", str(video), "-map", "0:v", "-c:v", "libx264", "-crf", "28",
                        "-preset", "slow", "-tune", "animation", "-pix_fmt", "yuv420p"])


def fingerprint(scene, flags):
    """What a chapter's picture and sound are made from, as far as this repository decides it: its
    own code, the code all chapters share (scenes.py without the other chapters), the kit, the
    model, the logo, and how it is rendered (the Manim version and the render flags). System
    tools such as LaTeX, Cairo, fonts and FFmpeg are assumed to stay put."""
    import manim
    module = ast.parse((HERE / "scenes.py").read_text())
    module.body = [node for node in module.body if not (
        isinstance(node, ast.ClassDef) and node.name in CHAPTERS and node.name != scene)]
    inputs = [ast.unparse(module).encode(), (HERE / "kit.py").read_bytes(),
              (HERE.parent / "eigencube.py").read_bytes(),
              (HERE.parent / "img" / "logo.svg").read_bytes(),
              " ".join([manim.__version__, *flags]).encode()]
    return hashlib.sha256(b"\0".join(inputs)).hexdigest()


def check_narration_audible(scene, video, media):
    """Fails the build if any spoken line came out silent; captions alone would hide that."""
    for cue in cues(media, scene):
        stats = subprocess.run(["ffmpeg", "-v", "info", "-ss", str(cue["start"]),
                                "-t", str(cue["end"] - cue["start"]), "-i", str(video),
                                "-af", "volumedetect", "-f", "null", "-"],
                               capture_output=True, text=True, check=True).stderr
        match = re.search(r"max_volume: (-?[\d.]+) dB", stats)
        peak = float(match.group(1)) if match else float("-inf")  # No audio track at all.
        if peak < -40:
            sys.exit(f"{scene}: narration is silent ({peak} dB) at {cue['start']:.1f}s: "
                     f"{cue['text']!r}")


def check_layout(scene, media):
    """Fails the build if, at the end of any animation, text overlapped what it must not (see
    Narrated.check_labels)."""
    collisions = json.loads((media / "cues" / f"{scene}.collisions.json").read_text())
    for time, name, other in collisions:
        print(f"{scene} at {time:.1f}s: {name!r} overlaps {other!r}", file=sys.stderr)
    if collisions:
        sys.exit(f"{scene}: {len(collisions)} text collision(s)")


def cues(media, scene):
    return json.loads((media / "cues" / f"{scene}.json").read_text())


def srt_time(seconds):
    millis = round(seconds * 1000)
    return "%02d:%02d:%02d,%03d" % (millis // 3600000, millis // 60000 % 60, millis // 1000 % 60,
                                    millis % 1000)


SPOKEN_NUMBER = r"(?:minus |plus )?(?:zero|one|two|three)"
SPOKEN_TUPLE = rf"(?i:\b{SPOKEN_NUMBER}(?:, (?:or |and )?{SPOKEN_NUMBER})+\b)"  # "one, zero, or one"


def caption_chunks(text):
    """Sentences, with long ones split at clause boundaries and short ones joined to a neighbor,
    so that each caption is comfortable to read. A spoken tuple ("one, zero, zero") stays on one
    caption."""
    # Its spaces become no-break spaces, which neither clause nor word breaks split at.
    text = re.sub(SPOKEN_TUPLE, lambda tuple_: tuple_[0].replace(" ", "\u00a0"), text)
    pieces = []
    for sentence in re.findall(r"[^.!?]+[.!?]*", text):
        for clause in re.split(r"(?<=[,:;]) +", sentence.strip()):
            # Still too long: break between words, into lines of even length. Filling each line
            # greedily would strand the clause's last word or two on a line of their own.
            while len(clause) > CAPTION_WIDTH:
                target = len(clause) / -(-len(clause) // CAPTION_WIDTH)
                spaces = [i for i, ch in enumerate(clause[:CAPTION_WIDTH + 1]) if ch == " "]
                if not spaces:
                    sys.exit(f"no place to break this caption: {clause!r}")
                cut = min(spaces, key=lambda i: abs(i - target))
                pieces.append(clause[:cut])
                clause = clause[cut + 1:]
            pieces.append(clause)
    chunks = [""]
    for piece in pieces:
        joined = f"{chunks[-1]} {piece}".strip()
        if len(joined) <= CAPTION_WIDTH and (len(chunks[-1]) < MIN_CAPTION or not
                                             re.search(r"[.!?]$", chunks[-1])):
            chunks[-1] = joined
        else:
            chunks.append(piece)
    # A short tail joins its predecessor; if the two don't fit on one line, they share it evenly,
    # breaking after punctuation where there is any, so as not to split a phrase.
    if len(chunks) > 1 and len(chunks[-1]) < MIN_CAPTION:
        joined = f"{chunks[-2]} {chunks[-1]}"
        if len(joined) <= CAPTION_WIDTH:
            chunks[-2:] = [joined]
        else:
            spaces = [i for i, ch in enumerate(joined) if ch == " "]
            breaks = [i for i in spaces if joined[i - 1] in ",;:.!?"] or spaces
            middle = min((i for i in breaks if max(i, len(joined) - i - 1) <= CAPTION_WIDTH),
                         key=lambda i: abs(i - len(joined) / 2), default=None)
            if middle is not None:
                chunks[-2:] = [joined[:middle], joined[middle + 1:]]
    return [chunk.replace("\u00a0", " ") for chunk in chunks]


def timed_chunks(cue):
    """Each caption chunk of a spoken line, timed by the words the voice actually spoke there.

    The voice speaks `speakable(text)`. Chunks and spoken words are both located in that text
    with a single cursor moving forward, so every word lands in the chunk that contains it."""
    spoken, bounds, position = cue["spoken"], [], 0
    chunks = caption_chunks(cue["text"])
    for chunk in chunks:
        position = spoken.index(speakable(chunk), position)
        bounds.append(position)
        position += len(speakable(chunk))
    starts, position = [None] * len(chunks), 0
    for word_start, _, word in cue["words"]:
        position = spoken.find(word, position)
        if position < 0:
            sys.exit(f"the voice reported {word!r}, which is not in the line: {cue['text']!r}")
        chunk = max(i for i, bound in enumerate(bounds) if bound <= position)
        if starts[chunk] is None:
            starts[chunk] = word_start
        position += len(word)
    for i in range(len(starts)):  # A chunk without a recognized word shares its neighbor's start.
        starts[i] = starts[i] if starts[i] is not None else (starts[i - 1] if i else cue["start"])
    ends = starts[1:] + [cue["end"]]  # A caption stays up until the next one starts.
    return list(zip(starts, ends, chunks))


def check_captions(entries):
    """Fails the build on captions that are out of order, too brief to read, or too long for one
    line."""
    for (a, b, text), (next_a, _, _) in zip(entries, entries[1:] + [(float("inf"), 0, "")]):
        if b - a < MIN_CAPTION_SECONDS or next_a < a or len(text) > CAPTION_WIDTH:
            sys.exit(f"caption at {a:.2f}s lasts {b - a:.2f}s, is out of order or is longer than "
                     f"{CAPTION_WIDTH} characters: {text!r}")


def captions(videos, media):
    entries, offset = [], 0.0
    for scene, video in zip(CHAPTERS, videos):
        for cue in cues(media, scene):
            entries += [(offset + a, offset + b, chunk) for a, b, chunk in timed_chunks(cue)]
        offset += duration(video)
    check_captions(entries)
    return "\n".join(f"{i}\n{srt_time(a)} --> {srt_time(b)}\n{text}\n"
                     for i, (a, b, text) in enumerate(entries, 1))


def check_av_lengths(film):
    """Fails the build unless the audio runs exactly as long as the picture."""
    lengths = dict(line.split(",") for line in subprocess.check_output(
        ["ffprobe", "-v", "error", "-show_entries", "stream=codec_type,duration", "-of", "csv=p=0",
         str(film)], text=True).split() if line.startswith(("video", "audio")))
    if abs(float(lengths["video"]) - float(lengths["audio"])) > 0.1:
        sys.exit(f"audio ({lengths['audio']}s) and video ({lengths['video']}s) differ in length")


def ffmpeg_to(path, args):
    """Runs ffmpeg into `path`, which only appears once ffmpeg has finished: an interrupted run
    could otherwise leave a short but valid file that later builds reuse as finished."""
    partial = path.with_name(f"{path.stem}.partial{path.suffix}")
    subprocess.run(["ffmpeg", "-v", "error", "-y", *args, str(partial)], check=True)
    partial.replace(path)


def concat_list(media, name, paths):
    """A listing of `paths` for ffmpeg's concat demuxer; returns its path."""
    listing = media / name
    quoted = (str(p).replace("'", "'\\''") for p in paths)  # The demuxer's quoting rule.
    listing.write_text("".join(f"file '{p}'\n" for p in quoted))
    return str(listing)


def contact_sheet(film, videos, media, tiles=36):
    """A grid of frames across the film, each taken as a spoken line ends: by then, the
    animations it narrates have settled, rather than being caught half-drawn."""
    ends, offset = [], 0.0
    for scene, video in zip(CHAPTERS, videos):
        ends += [offset + cue["end"] - 0.1 for cue in cues(media, scene)]
        offset += duration(video)
    frames = {round(ends[round(i * (len(ends) - 1) / (tiles - 1))] * FINAL_FPS) for i in range(tiles)}
    subprocess.run(["ffmpeg", "-v", "error", "-y", "-i", str(film), "-vf",
                    "select='%s',scale=384:216,tile=6x6" % "+".join(f"eq(n,{n})" for n in frames),
                    "-fps_mode", "passthrough", "-frames:v", "1", "-q:v", "4",
                    str(HERE / "contact-sheet.jpg")], check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--draft", action="store_true")
    args = parser.parse_args()
    flags, quality = (["-ql"], "480p15") if args.draft else (
        ["--resolution", "1920,1080", "--frame_rate", str(FINAL_FPS)], f"1080p{FINAL_FPS}")
    media = BUILD / ("draft" if args.draft else "final")
    for stage in ("fingerprints", "encoded"):
        (media / stage).mkdir(parents=True, exist_ok=True)

    def video(scene):
        """Renders the chapter if its inputs changed, checks it, and encodes its picture."""
        chapter = media / "videos" / "scenes" / quality / f"{scene}.mp4"
        stamp, current = media / "fingerprints" / scene, fingerprint(scene, flags)
        encoded = media / "encoded" / f"{scene}.mp4"
        if chapter.exists() and stamp.exists() and stamp.read_text() == current:
            print(f"{scene} is unchanged; reusing its earlier render", file=sys.stderr)
        else:
            stamp.unlink(missing_ok=True)  # Until it has rendered again, it is stale.
            render(scene, flags, media)
            encode(chapter, encoded)
            stamp.write_text(current)
        if not encoded.exists():  # A missing cache file is simply made again.
            encode(chapter, encoded)
        # Checked here, whether rendered just now or earlier: a chapter that failed its checks
        # stays on disk, and must not slip into a later build.
        check_narration_audible(scene, chapter, media)
        check_layout(scene, media)
        return chapter, encoded

    with ThreadPoolExecutor(min(len(CHAPTERS), os.cpu_count() or 1)) as pool:
        videos, pictures = zip(*pool.map(video, CHAPTERS))

    film = BUILD / ("draft.mp4" if args.draft else "eigencube-explainer.mp4")
    subtitles = film.with_suffix(".srt")
    subtitles.write_text(captions(videos, media))
    picture = sum(duration(v) for v in videos)
    # The soundtrack is mixed from all chapters at once (its loudness is the whole film's), and
    # kept by what it is made from: a change to the picture alone reuses it.
    # Each audio packet's timing as well as its content: a sound can move without changing.
    sounds = [subprocess.run(["ffmpeg", "-v", "error", "-i", str(v), "-map", "0:a", "-c", "copy",
                              "-f", "framecrc", "-"], capture_output=True, text=True,
                             check=True).stdout for v in videos]
    soundtrack = media / "soundtracks" / (hashlib.sha1(
        "".join([*sounds, str(picture)]).encode()).hexdigest() + ".m4a")
    if not soundtrack.exists():
        soundtrack.parent.mkdir(exist_ok=True)
        ffmpeg_to(soundtrack, [
            "-f", "concat", "-safe", "0", "-i", concat_list(media, "sounds.txt", videos),
            "-map", "0:a",
            # A chapter's audio ends at its last sound, which leaves gaps between chapters.
            # Players that ignore such gaps would run ahead of the picture, so they are filled
            # with silence. Then normalized to the loudness usual for video online (-16 LUFS).
            "-af", f"aresample=async=1:first_pts=0,apad=whole_dur={picture},"
                   "loudnorm=I=-16:TP=-1.5:LRA=11,aresample=48000",
            "-c:a", "aac", "-b:a", "160k"])
    subprocess.run(["ffmpeg", "-v", "error", "-y",
                    "-f", "concat", "-safe", "0",
                    "-i", concat_list(media, "pictures.txt", pictures),
                    "-i", str(soundtrack), "-i", str(subtitles),
                    "-map", "0:v", "-map", "1:a", "-map", "2",
                    "-c:v", "copy", "-c:a", "copy",  # Both encoded already.
                    "-c:s", "mov_text", "-metadata:s:s:0", "language=eng",
                    "-movflags", "+faststart", str(film)], check=True)
    check_av_lengths(film)
    if not args.draft:  # The committed contact sheet tracks the latest final cut.
        contact_sheet(film, videos, media)
    print(f"{film}  ({duration(film) / 60:.1f} min)")


if __name__ == "__main__":
    main()
