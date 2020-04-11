import re
import inspect


class FormatError(RuntimeError):
    pass


class Field:
    def __init__(self, name, ftype, default=None, optional=False, regex=None, min=None, max=None, test=None):
        self.name = name
        self.type = ftype
        self.default = default
        self.optional = optional
        if ftype != str and regex != None:
            raise FormatError("Cannot define regular expression constraint for non string field {0}".format(name))
        if type(regex) == str:
            regex = re.compile(regex)
        self.re = regex
        if self.type not in (int, float) and min is not None or max is not None:
            raise FormatError("Cannot define min or max constraint for non numeric field {0}".format(name))
        self.min = self.type(min) if min is not None else None
        self.max = self.type(max) if max is not None else None
        if test is not None and not callable(test):
            raise FormatError("Field {0} test function is not callable".format(name))

    def _getValueForType(self, name, ftype, source, context):
        if type(source) == ftype:
            value = source
        elif callable(ftype):
            try:
                if "context" in inspect.signature(ftype).parameters:
                    value = ftype(source, context=context)
                else:
                    value = ftype(source)
            except Exception as ex:
                raise ValueError("Cannot evaluate {0}: {1}".format(name, ex))
        elif type(ftype) == dict:
            if type(source) != str or source not in ftype:
                raise ValueError("Invalid value {0} for {1}".format(source, self.name))
            return ftype[source]
        if self.re is not None and not self.re.match(value):
            raise ValueError("Invalid value {0} for field {1}".format(value, self.name))
        if self.min is not None and value < self.min:
            raise ValueError("Invalid value {0} for field {1}".format(value, self.name))
        if self.max is not None and value > self.max:
            raise ValueError("Invalid value {0} for field {1}".format(value, self.name))
        return value

    def getValue(self, source, context=None):
        if source is None:
            if self.default is not None:
                return self.default
            return None

        ftype = self.type
        fname = self.name
        if type(ftype) != list:
            return self._getValueForType(fname, ftype, source, context)
        if len(ftype) > 1:
            return DictObject(ftype, source, context=context)
        ftype = ftype[0]
        if source is None:
            source = []
        if type(source) != list:
            raise ValueError("Field {0} must be a list".format(fname))
        value = []
        for nitem, item in enumerate(source):
            name = fname + "[{0}]".format(nitem)
            value.append(self._getValueForType(name, ftype, item, context))
        return value


class DictObject:
    def __init__(self, fields, value, context):
        if type(value) != dict:
            raise RuntimeError("Expected a dictionary - got {0}".format(type(value).__name__))
        used = set()
        for field in fields:
            if not isinstance(field, Field):
                raise FormatError("Invalid type {0} for field definition".format(type(field).__name__))
            fname = field.name
            if fname in used:
                raise FormatError("Field name {0} duplicated")
            if fname not in value and not field.optional:
                raise ValueError("Field {0} missing".format(fname))
            used.add(fname)
            fvalue = value.get(fname)
            fvalue = field.getValue(fvalue, context=context)
            setattr(self, fname, fvalue)
        unused = set(value) - used
        if len(unused) > 0:
            unusedfields = ",".join((str(x) for x in unused))
            raise ValueError("{0} contains unexpected fields {1}".format(type(self).__name__, unusedfields))
