import {API, getIdToken} from './actions';
import {getPassingFilterIndices} from './dataset_filter';
import {cacheValues, computeDerivedStats} from './VectorUtil';

function reshapeDistributionResult(distribution) {
  const results = [];
  distribution.forEach((distributionResult) => {
    const categories = distributionResult.categories;
    const dimension = distributionResult.name;
    for (let i = 0; i < distributionResult.values.length; i++) {
      results.push({
        dimension: dimension,
        name: categories[i],
        feature: distributionResult.values[i].name,
        mean: distributionResult.values[i].mean,
        percentExpressed: distributionResult.values[i].percentExpressed,
      });
    }
  });
  return results;
}

// convert sparse to dense and expand categorical values
function convertSparseAndCategoricalArrays(resultObject, ndim) {
  for (let key in resultObject) {
    const data = resultObject[key];
    if (data.indices) {
      // sparse
      const values = new Float32Array(ndim);
      for (let i = 0, n = data.indices.length; i < n; i++) {
        values[data.indices[i]] = data.values[i];
      }
      resultObject[key] = values;
    } else if (data.categories) {
      // category
      const values = new Array(ndim);
      for (let i = 0, n = data.values.length; i < n; i++) {
        values[i] = data.categories[data.values[i]];
      }
      resultObject[key] = values;
    }
  }
}

export class RestDataset {
  /**
   *
   * @param id Dataset id
   * param url Dataset URL
   * @param local Whether stats should be computed locally or by the server
   * @returns {Promise<void>}
   */
  init(id, url, local = true) {
    this.id = id;
    this.local = local;
    return Promise.resolve();
  }

  getSelectedIdsPromise(q, cachedData) {
    q.id = this.id;
    if (this.local) {
      const p =
        cachedData['index'] == null
          ? this.getDataPromise({values: {dimensions: ['index']}}, cachedData)
          : Promise.resolve();
      return p.then(() => {
        const selectedIndices = Array.from(
          getPassingFilterIndices(cachedData, q.filter),
        );
        const ids = [];
        const index = cachedData['index'];
        for (let i = 0, n = selectedIndices.length; i < n; i++) {
          ids.push(index[selectedIndices[i]]);
        }
        return {ids: ids};
      });
    } else {
      return fetch(API + '/selected_ids', {
        body: JSON.stringify(q),
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          Authorization: 'Bearer ' + getIdToken(),
        },
      }).then((result) => result.json());
    }
  }

  getFileUrl(file) {
    return (
      API +
      '/file?id=' +
      this.id +
      '&file=' +
      file +
      '&access_token=' +
      getIdToken()
    );
  }

  getJob(id, download = false) {
    return fetch(
      API +
        '/job?c=result&id=' +
        id +
        '&ds=' +
        this.id +
        (download ? '&dl=1' : ''),
      {
        headers: {Authorization: 'Bearer ' + getIdToken()},
      },
    ).then((response) => {
      if (download) {
        const disp = response.headers.get('content-disposition');
        let name = disp.substring(
          disp.indexOf('filename=') + 'filename='.length,
        );
        name = name.replaceAll('"', '');
        response.blob().then((val) => window.saveAs(val, name));
      } else {
        return response.json();
      }
    });
  }

  getJobParams(id) {
    return fetch(API + '/job?c=params&id=' + id, {
      headers: {Authorization: 'Bearer ' + getIdToken()},
    }).then((response) => {
      return response.json();
    });
  }

  getJobStatus(id) {
    return fetch(API + '/job?c=status&id=' + id, {
      headers: {Authorization: 'Bearer ' + getIdToken()},
    }).then((response) => {
      if (response.status === 404) {
        return null;
      }
      return response.json();
    });
  }

  getJobs(id) {
    return fetch(API + '/jobs?id=' + id, {
      headers: {Authorization: 'Bearer ' + getIdToken()},
    }).then((response) => {
      return response.json();
    });
  }

  deleteJob(id) {
    return fetch(API + '/job', {
      body: JSON.stringify({id: id}),
      method: 'DELETE',
      headers: {Authorization: 'Bearer ' + getIdToken()},
    });
  }

  getDataPromise(data, cachedData) {
    data.id = this.id;
    let dataSend = data;
    const local = this.local;

    if (local) {
      dataSend = {};
      //  ['stats', 'groupedStats', 'embedding', 'selection', 'values'];
      if (data.embedding || data.values) {
        dataSend.id = this.id;
        dataSend.embedding = data.embedding;
        dataSend.values = data.values;
      }
    }
    let jsonData = JSON.stringify(dataSend);
    let p =
      jsonData !== '{}'
        ? fetch(API + '/data', {
            body: jsonData,
            method: 'POST',
            headers: {
              'Content-Type': 'application/json',
              Authorization: 'Bearer ' + getIdToken(),
            },
          })
            .then((r) => r.json())
            .then((result) => {
              if (result.values) {
                convertSparseAndCategoricalArrays(
                  result.values,
                  this.schema.shape[0],
                );
              }
              if (result.layers) {
                for (const layer in result.layers) {
                  convertSparseAndCategoricalArrays(
                    result.layers[layer],
                    this.schema.shape[0],
                  );
                }
              }
              cacheValues(result, cachedData);
              return result;
            })
        : Promise.resolve({});
    return p.then((result) => {
      if (local) {
        computeDerivedStats(result, data, cachedData);
      } else {
        if (result.distribution) {
          result.distribution = reshapeDistributionResult(result.distribution);
        }
        if (result.selection && result.selection.distribution) {
          result.selection.distribution = reshapeDistributionResult(
            result.selection.distribution,
          );
        }
      }
      return result;
    });
  }

  getSchemaPromise() {
    if (this.schema != null) {
      return Promise.resolve(this.schema);
    }
    return fetch(API + '/schema?id=' + this.id, {
      headers: {Authorization: 'Bearer ' + getIdToken()},
    })
      .then((response) => {
        return response.json();
      })
      .then((result) => {
        this.schema = result;
        return result;
      });
  }
}
