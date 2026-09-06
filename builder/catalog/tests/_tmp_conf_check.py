def test_show(request):
    print("\ninifile:", request.config.inipath)
    print("rootdir:", request.config.rootpath)
    print("timeout:", repr(request.config.getini("timeout")))
    print("method :", repr(request.config.getini("timeout_method")))
