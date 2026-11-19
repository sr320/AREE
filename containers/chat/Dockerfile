# Ask AREE chat service (src/chat; docs/ask_aree.md#deploy-to-cloud-run).
#
# Build from the repository root:
#   docker build -f containers/chat/Dockerfile -t aree-chat .
# Run (the service reads its settings from AREE_CHAT_* variables):
#   docker run -p 8080:8080 -e ANTHROPIC_API_KEY=... aree-chat
#
# The image holds the code, the study registry and the published Explore-page tables,
# fetched at build time, so the chat answers from exactly the build the site serves.
# It does not hold the pipeline's reports/ or data/.
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# Locked versions (containers/chat/requirements.lock, from `uv pip compile pyproject.toml
# --extra chat`); the package itself is then installed without resolving anything.
COPY containers/chat/requirements.lock /tmp/requirements.lock
RUN pip install -r /tmp/requirements.lock

COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY registry ./registry
# Editable, so the code finds registry/ beside src/ as it does in a checkout.
RUN pip install --no-deps -e .

ARG SITE_DATA_URL=https://sr320.github.io/AREE/explore/data/
RUN aree fetch-site-data --dest /app/site-data --url "$SITE_DATA_URL"

RUN useradd --create-home --uid 10001 aree && chown -R aree /app
USER aree

ENV AREE_DATA_DIR=/app/site-data \
    AREE_CHAT_TRUST_PROXY=true \
    PORT=8080
EXPOSE 8080
CMD ["sh", "-c", "exec uvicorn --factory chat.server:app_from_env --host 0.0.0.0 --port \"$PORT\""]
