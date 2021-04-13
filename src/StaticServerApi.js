export class StaticServerApi {

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
