import requests
import json

query = {
  "size": 3000,
  "_source": ["name", "openAccess"],
  "query": {"match_all": {}}
}
r = requests.post("http://www.hipsci.org/lines/api/cellLine/_search", data = json.dumps(query))
dict = json.loads(r.text)

#Extract line id and openAccess status from the dict
lines = dict[u'hits'][u'hits']

for line in lines:
	print line[u'_id'] + "\t" + str(line[u'_source'][u'openAccess'])
