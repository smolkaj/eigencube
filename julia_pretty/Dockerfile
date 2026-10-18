FROM julia:1.13

RUN useradd --create-home --uid 10001 --user-group appuser \
    && install -d -o appuser -g appuser /app

WORKDIR /app

COPY --chown=appuser:appuser . .

USER appuser

# Install and precompile the port and its test environment, so the tests run
# without downloading or compiling packages.
RUN julia --project=. -e 'using Pkg; Pkg.instantiate()' \
    && julia --project=test -e 'using Pkg; Pkg.instantiate()'

HEALTHCHECK NONE

CMD ["julia", "--project=.", "-e", "using Pkg; Pkg.test()"]
