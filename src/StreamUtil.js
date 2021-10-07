import {getIdToken} from './actions';

function copyString(s) {
    return (' ' + s).substr(1);
}

export function streamingParseLines(url, readLine, done) {
    return new Promise((resolve, reject) => {
        fetch(url, {
            headers: {
                'Authorization': 'Bearer ' + getIdToken(),
                'Content-Type': 'application/json'
            }
        }).then(function (response) {
            if (response.ok) {
                const reader = response.body.getReader();
                const textDecoder = new TextDecoder();
                let skipLF = false;
                let text = '';
                reader.read().then(function processResult(result) {
                    // result contains a value which is an array of Uint8Array
                    text += (result.done ? '' : textDecoder(result.value, 0, result.value.length));
                    let start = 0;
                    // TODO no need to search previous chunk of text
                    for (let i = 0, length = text.length; i < length; i++) {
                        let c = text[i];
                        if (skipLF && c === '\n') {
                            start++;
                            skipLF = false;
                        } else if (c === '\n' || c === '\r') {
                            skipLF = c === '\r'; // \r\n windows line ending
                            const s = copyString(text.substring(start, i));
                            readLine(s);
                            start = i + 1;
                        } else {
                            skipLF = false;
                        }
                    }
                    text = start < text.length ? text.substring(start) : '';
                    if (!result.done) {
                        return reader.read().then(processResult);
                    } else {
                        if (text !== '' && text !== '\r') {
                            readLine(text);
                        }
                        resolve(done());
                    }
                });
            } else {
                reject('Network error');
            }
        }).catch(function (error) {
            reject(error);
        });
    });


}