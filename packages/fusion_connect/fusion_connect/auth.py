"""OAuth2 authentication against Autodesk Platform Services (APS).

Supports the 2-legged (client credentials) flow, which is the right choice for
server-to-server access to your own Autodesk account's data. For 3-legged user
flows you can pass a pre-obtained access token to APSAuth.from_token instead.
"""

import os
import time
import threading

import requests


APS_TOKEN_URL = 'https://developer.api.autodesk.com/authentication/v2/token'

DEFAULT_SCOPES = (
    'data:read',
    'data:write',
    'data:create',
    'bucket:read',
    'bucket:create',
)


class AuthError(RuntimeError):
    pass


class APSAuth:
    """Manages an APS access token, refreshing it before expiry."""

    def __init__(self, client_id, client_secret, scopes=DEFAULT_SCOPES,
                 token_url=APS_TOKEN_URL, refresh_skew=60):
        if not client_id or not client_secret:
            raise AuthError('client_id and client_secret are required')
        self._client_id = client_id
        self._client_secret = client_secret
        self._scopes = tuple(scopes)
        self._token_url = token_url
        self._refresh_skew = refresh_skew
        self._lock = threading.Lock()
        self._access_token = None
        self._expires_at = 0.0

    @classmethod
    def from_env(cls, scopes=DEFAULT_SCOPES):
        # APS_CLIENT_ID / APS_CLIENT_SECRET are the canonical env var names
        # documented by Autodesk for APS apps.
        client_id = os.environ.get('APS_CLIENT_ID') or os.environ.get('FORGE_CLIENT_ID')
        client_secret = os.environ.get('APS_CLIENT_SECRET') or os.environ.get('FORGE_CLIENT_SECRET')
        return cls(client_id, client_secret, scopes=scopes)

    @classmethod
    def from_token(cls, access_token, expires_in=3600):
        """Wrap a token obtained out-of-band (e.g. a 3-legged user token)."""
        obj = cls.__new__(cls)
        obj._client_id = None
        obj._client_secret = None
        obj._scopes = ()
        obj._token_url = APS_TOKEN_URL
        obj._refresh_skew = 60
        obj._lock = threading.Lock()
        obj._access_token = access_token
        obj._expires_at = time.time() + float(expires_in)
        return obj

    def token(self):
        with self._lock:
            if self._access_token and time.time() < self._expires_at - self._refresh_skew:
                return self._access_token
            if not self._client_id:
                # Token was supplied via from_token and we cannot refresh it.
                if self._access_token:
                    return self._access_token
                raise AuthError('access token has expired and no client credentials are available to refresh it')
            self._fetch()
            return self._access_token

    def _fetch(self):
        response = requests.post(
            self._token_url,
            data={
                'grant_type': 'client_credentials',
                'scope': ' '.join(self._scopes),
            },
            auth=(self._client_id, self._client_secret),
            headers={'Accept': 'application/json'},
            timeout=30,
        )
        if response.status_code != 200:
            raise AuthError(
                'APS token request failed ({}): {}'.format(response.status_code, response.text)
            )
        payload = response.json()
        self._access_token = payload['access_token']
        self._expires_at = time.time() + float(payload.get('expires_in', 3600))
