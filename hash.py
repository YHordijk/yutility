import hashlib
import json
import jsonpickle


def hash(obj):
    '''
    Hash any python object using SHA-256. If the object is a dictionary use JSON to
    put it in a standard format. If there are non-picklable objects
    use jsonpickle to pickle them anyway.
    '''
    if isinstance(obj, dict):
        try:
            obj = json.dumps(obj, sort_keys=True, indent=4)
        except:
            obj = jsonpickle.encode(obj)

    return hashlib.sha256(obj.encode('utf-8')).hexdigest()
