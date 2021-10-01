FROM node:14.17.1-buster-slim
SHELL ["/bin/bash", "-c"]
RUN apt-get update && \
    apt-get install --no-install-recommends -y python3-dev python3-pip python3-wheel git && \
    ln -s /usr/bin/python3 /usr/bin/python
COPY . /cirrocumulus
WORKDIR /cirrocumulus
RUN yarn global add typescript
RUN yarn install
RUN yarn build
RUN python -m pip install --upgrade pip
RUN python -m pip install setuptools
RUN python -m pip install .

EXPOSE 3000
CMD ["cirro", "serve", "--db_uri", "", "--bind", "0.0.0.0:3000"]
