const authScopes = ['email'
    // 'profile',
    // 'https://www.googleapis.com/auth/userinfo.profile',
    // 'https://www.googleapis.com/auth/contacts.readonly',
    // 'https://www.googleapis.com/auth/devstorage.full_control',
];

export function GoggleAuth() {

    this.signOut = function () {
        return window.gapi.auth2.getAuthInstance().signOut();
    };

    this.signIn = function () {
        return window.gapi.auth2.getAuthInstance().signIn();
    };

    this.getEmail = function () {
        const isSignedIn = window.gapi.auth2.getAuthInstance().isSignedIn.get();
        if (isSignedIn) {
            return window.gapi.auth2.getAuthInstance().currentUser.get().getBasicProfile().getEmail();
        }
    };

    // this.getAccessToken = function () {
    //     return window.gapi.auth2.getAuthInstance().currentUser.get().getAuthResponse().access_token;
    // };

    this.getIdToken = function () {
        return typeof window.gapi !== 'undefined' ? window.gapi.auth2.getAuthInstance().currentUser.get().getAuthResponse().id_token : '';
    };

    this.init = function (serverInfo) {

        return new Promise((resolve) => {
            const script = document.createElement('script');
            script.type = 'text/javascript';
            script.src = 'https://apis.google.com/js/api.js';

            script.onload = (e) => {
                window.gapi.load('client:auth2', () => {
                    window.gapi.client.init({
                        clientId: serverInfo.clientId, scope: authScopes.join(' ')
                    }).then(() => {
                        resolve();
                    });
                });
            };
            document.getElementsByTagName('head')[0].appendChild(script);
        });
    };
}
