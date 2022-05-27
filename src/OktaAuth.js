import {OktaAuth as OktaAuthClient} from '@okta/okta-auth-js';

export function OktaAuth() {
  const _this = this;
  _this.idToken = null;
  this.signOut = function () {
    return _this.authClient.signOut().then(() => {
      _this.idToken = null;
    });
  };

  this.signIn = function () {
    return _this.authClient.token
      .getWithPopup({
        nonce: null,
        scopes: ['email', 'openid', 'profile'],
      })
      .then((res) => {
        _this.idToken = res.tokens.idToken;
        _this.authClient.tokenManager.add('idToken', res.tokens.idToken);
        return res;
      })
      .catch((err) => {
        console.log(err);
      });
  };

  this.getEmail = function () {
    return _this.idToken ? _this.idToken.claims.email : null;
  };

  this.hasExpired = function () {
    return _this.authClient.tokenManager.hasExpired(_this.idToken);
  };

  this.getIdToken = function () {
    return _this.idToken.idToken;
  };

  this.init = function (authInfo) {
    const config = {
      issuer: authInfo.issuer,
      clientId: authInfo.clientId,
      postLogoutRedirectUri:
        window.location.protocol + '//' + window.location.host,
      pkce: true,
    };

    return new Promise((resolve) => {
      _this.authClient = new OktaAuthClient(config);
      const handler = (authState) => {
        _this.idToken = authState.idToken;
      };
      _this.authClient.authStateManager.subscribe(handler);
      _this.authClient.start(); // start the service
      if (_this.authClient.isLoginRedirect()) {
        _this.authClient.token.parseFromUrl().then((data) => {
          const {idToken} = data.tokens;
          _this.idToken = idToken;
          resolve();
        });
      } else {
        _this.authClient.tokenManager
          .get('idToken')
          .then((idToken) => {
            if (idToken) {
              _this.idToken = idToken;
              return _this.authClient.tokenManager.hasExpired(idToken);
            }
            return true;
          })
          .then((expired) => {
            if (expired) {
              _this.idToken = null;
            }
            resolve();
          });
      }
    });
  };
}
