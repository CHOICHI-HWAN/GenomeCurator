# Virus Genome Metadata Curator - same conda environment as environment.yml, multi-arch (linux/amd64, linux/arm64)
FROM condaforge/miniforge3:latest
WORKDIR /app
COPY environment.yml pyproject.toml README.md LICENSE ./
COPY genomecurator ./genomecurator
# build the environment from environment.yml (its last line installs this package with `pip -e .`)
RUN mamba env create -f environment.yml && mamba clean -afy
ENV GC_DATA_DIR=/data PATH=/opt/conda/envs/genomecurator/bin:$PATH
VOLUME /data
EXPOSE 8001
# the app has no login: publish the port to localhost only (docker run -p 127.0.0.1:8001:8001 ...)
CMD ["python", "-m", "genomecurator", "serve", "--host", "0.0.0.0", "--port", "8001"]
