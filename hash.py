import hashlib
import json
import jsonpickle


def hash(obj, use_json=True):
    if use_json:
        try:
            obj = json.dumps(obj, sort_keys=True, indent=4)
        except:
            obj = jsonpickle.encode(obj)

    return hashlib.sha256(obj.encode('utf-8')).hexdigest()
