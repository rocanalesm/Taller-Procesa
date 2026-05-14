# fusion_connect

A small Python client for the Autodesk Platform Services (APS, formerly Forge)
REST API used by Autodesk Fusion.  It handles OAuth2 token management and
exposes the most common Data Management endpoints (hubs, projects, folders,
items, versions) so notebooks and scripts can read Fusion files from your
account.

> **Scope** — this package is read-oriented and intentionally tiny.  It is not
> the in-process Fusion 360 Python API (which only runs inside the desktop
> Fusion application).  If you need to script Fusion locally, use the
> built-in [Fusion 360 API](https://help.autodesk.com/view/fusion360/ENU/?guid=GUID-A92A4B10-3781-4925-94C6-47DA85A4F65A)
> instead.

## Install

From the repository root:

```bash
cd packages/fusion_connect
pip install .
```

## Get APS credentials

1. Go to https://aps.autodesk.com/ and create an application.
2. Copy the **Client ID** and **Client Secret**.
3. Export them as environment variables:

```bash
export APS_CLIENT_ID="..."
export APS_CLIENT_SECRET="..."
```

For account data (Fusion Team hub, Personal hub) the 2-legged client-credentials
flow used by `fusion_connect` requires that your APS app be **provisioned** for
the hub.  See Autodesk's docs on provisioning an APS app to a Fusion Team / BIM
360 / ACC hub.

## Quick start

```python
from fusion_connect import APSAuth, FusionClient

auth = APSAuth.from_env()              # reads APS_CLIENT_ID / APS_CLIENT_SECRET
fusion = FusionClient(auth)

for hub in fusion.list_hubs():
    print(hub["id"], hub["attributes"]["name"])
    for project in fusion.list_projects(hub["id"]):
        print("  ", project["id"], project["attributes"]["name"])
```

Walk into a project and list top-level folders:

```python
hub_id = "b.xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
project_id = "a.yyyyyyyy"

for folder in fusion.get_top_folders(hub_id, project_id):
    print(folder["id"], folder["attributes"]["displayName"])
```

Resolve the tip version of an item and download it:

```python
tip = fusion.get_tip_version(project_id, item_id)
version_id = tip["data"]["id"]
fusion.download_version(project_id, version_id, "/tmp/part.f3d")
```

## Using a 3-legged user token

If you already have a 3-legged access token (e.g. from a web app's user login),
wrap it directly:

```python
auth = APSAuth.from_token("<access-token>", expires_in=3600)
fusion = FusionClient(auth)
```

## API surface

| Method | Description |
| --- | --- |
| `list_hubs()` | All APS hubs visible to the credentials. |
| `list_projects(hub_id)` | Projects inside a hub. |
| `get_project(hub_id, project_id)` | Project metadata. |
| `get_top_folders(hub_id, project_id)` | Project root folders. |
| `list_folder_contents(project_id, folder_id, filter_type=None)` | Items + sub-folders in a folder. |
| `get_item(project_id, item_id)` | Item metadata. |
| `list_item_versions(project_id, item_id)` | Every version of an item. |
| `get_tip_version(project_id, item_id)` | Latest version of an item. |
| `download_version(project_id, version_id, dest_path)` | Stream a version's binary to disk. |
