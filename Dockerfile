FROM python:3.13-slim

LABEL com.azure.containerizationassist.createdby="containerization-assist"

ENV PYTHONUNBUFFERED=1 \
    SDL_VIDEODRIVER=dummy

WORKDIR /app

COPY requirements.txt .
RUN python -m pip install --no-cache-dir -r requirements.txt \
    && useradd --create-home --uid 10001 --user-group appuser

COPY --chown=appuser:appuser . .

HEALTHCHECK NONE

USER appuser

CMD ["python", "-m", "unittest", "discover", "tests"]
