# Test         Grid    PEs        Sets    BFB-compare
smoke          tx1m    32x1      run10day,gridc
smoke          tx1m    32x1      run10day
restart        tx1m    16x1      gridc  smoke_tx1m_32x1_gridc_run10day
restart        tx1m    16x1      debug,gridc
