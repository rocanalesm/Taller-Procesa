"""Thin wrapper over the Autodesk Platform Services (APS) Data Management API.

Fusion (Team / Personal) files live inside APS hubs and projects. The methods
below cover the read paths most users need: listing hubs and projects, walking
folders, and resolving item versions. Write/upload paths are deliberately out
of scope here.
"""

import urllib.parse

import requests


APS_API_BASE = 'https://developer.api.autodesk.com'


class APSError(RuntimeError):
    def __init__(self, status_code, message, payload=None):
        super().__init__('APS API error {}: {}'.format(status_code, message))
        self.status_code = status_code
        self.payload = payload


class FusionClient:
    """Read-oriented client for Fusion data hosted in APS."""

    def __init__(self, auth, base_url=APS_API_BASE, session=None, timeout=30):
        self._auth = auth
        self._base_url = base_url.rstrip('/')
        self._session = session or requests.Session()
        self._timeout = timeout

    def _request(self, method, path, **kwargs):
        url = path if path.startswith('http') else '{}{}'.format(self._base_url, path)
        headers = kwargs.pop('headers', {}) or {}
        headers.setdefault('Authorization', 'Bearer {}'.format(self._auth.token()))
        headers.setdefault('Accept', 'application/vnd.api+json')
        response = self._session.request(
            method, url, headers=headers, timeout=self._timeout, **kwargs
        )
        if not response.ok:
            try:
                payload = response.json()
            except ValueError:
                payload = response.text
            raise APSError(response.status_code, response.reason, payload)
        if response.status_code == 204 or not response.content:
            return None
        return response.json()

    def _paginated(self, path, params=None):
        """Follow JSON:API `links.next` pages and yield `data` items."""
        next_url = path
        next_params = params
        while next_url:
            payload = self._request('GET', next_url, params=next_params)
            for item in payload.get('data', []) or []:
                yield item
            next_url = (payload.get('links') or {}).get('next', {}).get('href')
            next_params = None  # The `next` link already encodes params.

    # ----- Hubs / projects -------------------------------------------------

    def list_hubs(self):
        """Return all hubs (Fusion Teams, A360 personal, BIM360) visible to the credentials."""
        return list(self._paginated('/project/v1/hubs'))

    def list_projects(self, hub_id):
        return list(self._paginated('/project/v1/hubs/{}/projects'.format(hub_id)))

    def get_project(self, hub_id, project_id):
        return self._request('GET', '/project/v1/hubs/{}/projects/{}'.format(hub_id, project_id))

    def get_top_folders(self, hub_id, project_id):
        payload = self._request(
            'GET',
            '/project/v1/hubs/{}/projects/{}/topFolders'.format(hub_id, project_id),
        )
        return payload.get('data', []) if payload else []

    # ----- Folders / items / versions --------------------------------------

    def list_folder_contents(self, project_id, folder_id, filter_type=None):
        params = {}
        if filter_type:
            params['filter[type]'] = filter_type
        return list(self._paginated(
            '/data/v1/projects/{}/folders/{}/contents'.format(project_id, folder_id),
            params=params or None,
        ))

    def get_item(self, project_id, item_id):
        return self._request('GET', '/data/v1/projects/{}/items/{}'.format(project_id, item_id))

    def list_item_versions(self, project_id, item_id):
        return list(self._paginated(
            '/data/v1/projects/{}/items/{}/versions'.format(project_id, item_id)
        ))

    def get_tip_version(self, project_id, item_id):
        return self._request(
            'GET',
            '/data/v1/projects/{}/items/{}/tip'.format(project_id, item_id),
        )

    def download_version(self, project_id, version_id, dest_path):
        """Download the storage payload behind a version to a local file."""
        version = self._request(
            'GET',
            '/data/v1/projects/{}/versions/{}'.format(
                project_id, urllib.parse.quote(version_id, safe='')
            ),
        )
        storage = (((version or {}).get('data') or {}).get('relationships') or {}).get('storage')
        if not storage or not storage.get('data'):
            raise APSError(404, 'version has no associated storage object')
        storage_id = storage['data']['id']  # e.g. urn:adsk.objects:os.object:<bucket>/<object>
        try:
            _, object_part = storage_id.split('urn:adsk.objects:os.object:', 1)
            bucket_key, object_key = object_part.split('/', 1)
        except ValueError:
            raise APSError(500, 'unrecognized storage URN: {}'.format(storage_id))

        signed = self._request(
            'GET',
            '/oss/v2/buckets/{}/objects/{}/signeds3download'.format(
                bucket_key, urllib.parse.quote(object_key, safe='')
            ),
            headers={'Accept': 'application/json'},
        )
        download_url = signed['url']
        with self._session.get(download_url, stream=True, timeout=self._timeout) as r:
            r.raise_for_status()
            with open(dest_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=1 << 16):
                    if chunk:
                        f.write(chunk)
        return dest_path
