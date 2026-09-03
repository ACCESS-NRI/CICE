# Test         Grid          PEs        Sets    BFB-compare
smoke          access_100km  32x1      run10day,gridc
smoke          access_100km  32x1      run10day
restart        access_100km  16x1      gridc  smoke_access_100km_32x1_gridc_run10day
restart        access_100km  16x1      debug,gridc