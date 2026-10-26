#!/usr/bin/env bash
# One-time Google Cloud setup for the Ask AREE chat service (docs/ask_aree.md#deploy-to-cloud-run).
#
# Run it in Google Cloud Shell (the terminal button in the Cloud console), which has
# gcloud and is already signed in:
#
#   bash setup_cloud_run.sh PROJECT_ID [REGION]
#
# It creates, and on a re-run leaves in place:
#   - the APIs Cloud Run, Artifact Registry, Secret Manager and IAM credentials;
#   - an Artifact Registry repository "aree" for the image;
#   - a Secret Manager secret "aree-anthropic-api-key" (you paste the key; it is not echoed);
#   - a runtime service account that may read only that secret;
#   - a deploy service account that GitHub Actions in sr320/AREE (and no other repository)
#     can act as through Workload Identity Federation, so no key file exists anywhere.
# Then it prints the GitHub repository variables the deploy workflow reads.
set -euo pipefail

PROJECT="${1:?usage: bash setup_cloud_run.sh PROJECT_ID [REGION]}"
REGION="${2:-us-west1}"
GITHUB_REPO="sr320/AREE"
SECRET="aree-anthropic-api-key"
RUNTIME_SA="aree-chat-runtime@${PROJECT}.iam.gserviceaccount.com"
DEPLOY_SA="aree-chat-deployer@${PROJECT}.iam.gserviceaccount.com"
POOL="github"
PROVIDER="aree"

exists() { "$@" >/dev/null 2>&1; }

gcloud config set project "$PROJECT" >/dev/null
PROJECT_NUMBER="$(gcloud projects describe "$PROJECT" --format='value(projectNumber)')"

echo "== APIs"
gcloud services enable run.googleapis.com artifactregistry.googleapis.com \
  secretmanager.googleapis.com iamcredentials.googleapis.com sts.googleapis.com

echo "== Artifact Registry repository 'aree' in $REGION"
exists gcloud artifacts repositories describe aree --location "$REGION" ||
  gcloud artifacts repositories create aree --repository-format=docker --location "$REGION" \
    --description="Ask AREE chat service images"

echo "== Secret '$SECRET'"
if exists gcloud secrets describe "$SECRET"; then
  echo "exists; to replace the key: gcloud secrets versions add $SECRET --data-file=-"
else
  read -rsp "Paste the Anthropic API key for the Ask AREE workspace (not shown): " KEY; echo
  [ -n "$KEY" ] || { echo "No key entered." >&2; exit 1; }
  printf %s "$KEY" | gcloud secrets create "$SECRET" --replication-policy=automatic --data-file=-
  unset KEY
fi

echo "== Service accounts"
exists gcloud iam service-accounts describe "$RUNTIME_SA" ||
  gcloud iam service-accounts create aree-chat-runtime --display-name="Ask AREE chat (runtime)"
exists gcloud iam service-accounts describe "$DEPLOY_SA" ||
  gcloud iam service-accounts create aree-chat-deployer --display-name="Ask AREE chat (GitHub deploy)"

# The runtime account reads the one secret; the deploy account deploys Cloud Run services,
# pushes images, and runs the service as the runtime account.
gcloud secrets add-iam-policy-binding "$SECRET" --quiet \
  --member="serviceAccount:$RUNTIME_SA" --role=roles/secretmanager.secretAccessor >/dev/null
for role in roles/run.admin roles/artifactregistry.writer; do
  gcloud projects add-iam-policy-binding "$PROJECT" --quiet --condition=None \
    --member="serviceAccount:$DEPLOY_SA" --role="$role" >/dev/null
done
gcloud iam service-accounts add-iam-policy-binding "$RUNTIME_SA" --quiet \
  --member="serviceAccount:$DEPLOY_SA" --role=roles/iam.serviceAccountUser >/dev/null

echo "== Workload Identity Federation for GitHub Actions in $GITHUB_REPO"
exists gcloud iam workload-identity-pools describe "$POOL" --location=global ||
  gcloud iam workload-identity-pools create "$POOL" --location=global --display-name="GitHub Actions"
exists gcloud iam workload-identity-pools providers describe "$PROVIDER" --location=global \
  --workload-identity-pool="$POOL" ||
  gcloud iam workload-identity-pools providers create-oidc "$PROVIDER" --location=global \
    --workload-identity-pool="$POOL" --display-name="sr320/AREE" \
    --issuer-uri="https://token.actions.githubusercontent.com" \
    --attribute-mapping="google.subject=assertion.sub,attribute.repository=assertion.repository" \
    --attribute-condition="assertion.repository == '$GITHUB_REPO'"
gcloud iam service-accounts add-iam-policy-binding "$DEPLOY_SA" --quiet \
  --role=roles/iam.workloadIdentityUser \
  --member="principalSet://iam.googleapis.com/projects/$PROJECT_NUMBER/locations/global/workloadIdentityPools/$POOL/attribute.repository/$GITHUB_REPO" \
  >/dev/null

WIF="projects/$PROJECT_NUMBER/locations/global/workloadIdentityPools/$POOL/providers/$PROVIDER"
cat <<EOF

Done. Set these GitHub repository variables (Settings > Secrets and variables > Actions >
Variables), or from a machine with the GitHub CLI:

gh variable set GCP_PROJECT_ID   --repo $GITHUB_REPO --body "$PROJECT"
gh variable set GCP_REGION       --repo $GITHUB_REPO --body "$REGION"
gh variable set GCP_WIF_PROVIDER --repo $GITHUB_REPO --body "$WIF"
gh variable set GCP_DEPLOY_SA    --repo $GITHUB_REPO --body "$DEPLOY_SA"
gh variable set GCP_RUNTIME_SA   --repo $GITHUB_REPO --body "$RUNTIME_SA"

None of these is a secret. Then run the "Deploy Ask AREE chat service" workflow
(Actions tab, or: gh workflow run chat-deploy.yml --repo $GITHUB_REPO).
EOF
