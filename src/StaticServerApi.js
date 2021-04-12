// general server stuff

export class StaticServerApi {


    isStaticApi() {
        return true;
    }

    getDatasetsPromise() {
        return fetch('/datasets.json')
            .then(response => {
                const results = response.json();
                for (let i = 0; i < results.length; i++) {
                    results[i].access = 'direct';
                }
                return results;
            });
    }


}
