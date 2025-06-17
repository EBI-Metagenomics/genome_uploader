class EnaQueryError(Exception):
    pass

class InvalidAccessionError(EnaQueryError):
    pass

class EnaEmptyResponseError(EnaQueryError):
    pass

class EnaParseError(EnaQueryError):
    pass

class NoDataException(EnaQueryError):
    pass