const authScopes = [
  'email', // 'profile',
  // 'https://www.googleapis.com/auth/userinfo.profile',
  // 'https://www.googleapis.com/auth/contacts.readonly',
  // 'https://www.googleapis.com/auth/devstorage.full_control',
];

const LOCAL_STORAGE_KEY = 'google-token-storage';

export function GoggleAuth(props) {
  let credential = null;
  let user = null;

  this.getIdToken = function () {
    return credential;
  };

  this.getEmail = function () {
    return credential && user ? user.email : null;
  };

  this.signIn = function () {
    // never gets called
  };

  this.signOut = function () {
    return new Promise((resolve) => {
      window.google.accounts.id.revoke(user.sub, (done) => {
        window.google.accounts.id.disableAutoSelect();
        credential = null;
        user = null;
        delete window.localStorage[LOCAL_STORAGE_KEY];
        resolve();
      });
    });
  };

  this.init = function (authInfo, api, callback) {
    return new Promise((resolve) => {
      function handleCredentialResponse(response) {
        credential = response.credential; //  ID token as a base64-encoded JSON Web Token (JWT) string.
        // get user from server
        api
          .getUserPromise()
          .then((_user) => {
            user = _user;
            window.localStorage[LOCAL_STORAGE_KEY] = credential;
            callback();
          })
          .catch((err) => {
            delete window.localStorage[LOCAL_STORAGE_KEY];
            console.log(err);
          });
      }

      function onGapiLoad() {
        window.google.accounts.id.initialize({
          client_id: authInfo.clientId, // itp_support: true, // for one tap
          scope: authScopes.join(' '),
          callback: handleCredentialResponse,
        });

        window.google.accounts.id.renderButton(
          document.getElementById('login-button'),
          {
            theme: 'outline',
            size: 'large',
          },
        );
        if (LOCAL_STORAGE_KEY in window.localStorage) {
          credential = window.localStorage[LOCAL_STORAGE_KEY];
          api
            .getUserPromise()
            .then((_user) => {
              user = _user;
            })
            .catch((err) => {
              console.log('Unable to get user');
              delete window.localStorage[LOCAL_STORAGE_KEY];
            })
            .finally(() => resolve());
        } else {
          resolve();
        }
      }

      function addScript(src, onLoad) {
        const script = document.createElement('script');
        script.type = 'text/javascript';
        script.async = true;
        script.defer = true;
        script.src = src;

        script.onload = (e) => {
          onLoad();
        };

        document.getElementsByTagName('head')[0].appendChild(script);
      }

      addScript('https://accounts.google.com/gsi/client', onGapiLoad);
    });
  };
}
