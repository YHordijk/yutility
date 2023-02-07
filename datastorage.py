from yutility.ytypes import check_hints, List


def Column:
	def __init__(self, key, data, unit=None, labels=None):
		self.key = key
		self.data = data
		self.unit = unit
		self.labels = labels or []


class DataStorage:
	def __init__(self, keys, data, units=None, labels=None):
		self.keys = keys
		self.data = data
		self.units = units or {}
		self.labels = labels or {}

	## KEY METHODS
    @check_hints
    def check_key(self, key: Key):
        return all([key in lst for lst in [self.values, self.units, self.labels, self.tags]])

    @check_hints
    def get_keys(self, tags: List(str) = None) -> List(str):
        return list(self.values.keys())

    @check_hints
    def get_subkeys(self, headkey: Key):
        if headkey == '*':
            return self.get_keys()

        headkey = self.parse_key(headkey)
        subkeys = []
        for key in self.get_keys():
            key = self.parse_key(key)
            for keypart, headkeypart in zip(key, headkey):
                if keypart != headkeypart:
                    break
            else:
                subkeys.append(self.key_to_string(key))
        return subkeys

    @check_hints
    def parse_key(self, key: Key) -> List(str):
        '''
        Parses a key from specific format.

        Keys are given as key1:key2:...:keyn or as an iterable [key1, key2, ..., keyn].
        This function ensures it is given in list format.
        '''
        if isinstance(key, str):
            return key.split(':')
        else:
            return key

    @check_hints
    def key_to_string(self, key: Key) -> str:
        '''
        Returns `key` in string format with key parts separated by `:`.
        '''
        key = self.parse_key(key)
        return ':'.join(key)


     ## RETRIEVAL
    @check_hints
    def __contains__(self, key: str) -> bool:
        '''
        Checks if a key is present.
        '''
        try:
            self.get(key)
            return True
        except KeyError:
            return False

    @check_hints
    def get_save(self, key: Key = '*', tags: List(str) = None, from_dict: dict = None) -> Either(Value, dict):
        '''
        Returns Datum object corresponding to a key.

        Parameters
        ----------
            key : Either(str, List(str))
                Key to access. Key are given as key1:key2:...:keyn.
                Key can be given as *, which will return all data.
            tags : List(str), optional
                List of tags that the Datum(s) must have.

        Returns
        -------
            Either(Datum, dict)
                Datum or a dictionary of Datums corresponding to the key.

        Raises
        ------
            KeyError 
                If the key could not be found.

        See Also
        --------
            ychem.results.datamanger.Datum
        '''
        if from_dict is None:
            from_dict = self.values

        subkeys = self.get_subkeys(key)

        if len(subkeys) == 0:
            raise KeyError(f'Could not find key {key}')

        if len(subkeys) == 1 and key == subkeys[0]:
            return from_dict[key]

        for subkey in subkeys:
            if not self.check_key(subkey):
                raise KeyError(f'Could not find key {subkey}')
        data_lst = [self.parse_key(subkey.replace(key + ':', '')) + [from_dict[subkey]] for subkey in subkeys]
        data = dictfunc.list_to_dict(data_lst)
        return data

        if tags is not None:
            return self.filter_data(data, tags)

        return data

    @check_hints
    def get(self, key: Key = '*', tags: List(str) = None) -> Either(Value, dict):
        try:
            return self.get_save(key=key, tags=tags)
        except:
            return None