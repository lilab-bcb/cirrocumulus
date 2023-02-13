set -e

mongosh <<EOF
db = db.getSiblingDB('$CIRRO_DATABASE')

db.createUser({
  user: '$CIRRO_DATABASE_USER',
  pwd: '$CIRRO_DATABASE_PASSWORD',
  roles: [{ role: 'readWrite', db: '$CIRRO_DATABASE' }],
});

EOF
