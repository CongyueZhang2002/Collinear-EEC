    hq0 = hq0_raw
    hg0 = hg0_raw
    hq1 = (hq1_raw - 2.3333333333333335*hg0_raw*LQ_ini + 5.555555555555555*hq0_raw*LQ_ini)
    hg1 = (hg1_raw + 11.733333333333334*hg0_raw*LQ_ini - 1.5555555555555556*hq0_raw*LQ_ini)
    hq2 = (
      hq2_raw + LQ_ini*(-2.333333333333333*hg1_raw + 32.43872728705069*hq0_raw + 13.222222222222221*hq1_raw + 
      hg0_raw*(19.78716621742935 - 29.114814814814814*LQ_ini) + 38.54320987654321*hq0_raw*LQ_ini)
    )
    hg2 = (
      hg2_raw + LQ_ini*(19.400000000000002*hg1_raw - 44.88570687115846*hq0_raw - 1.5555555555555556*hq1_raw - 19.40987654320988*hq0_raw*LQ_ini + 
      hg0_raw*(-79.94137331684473 + 115.62814814814814*LQ_ini))
    )
    dhq1 = (dhq1_raw + 0.6611111111111111*hg0_raw*LQ_ini + 1.6433520602275775*hq0_raw*LQ_ini)
    dhg1 = (dhg1_raw + 5.3425421355120495*hg0_raw*LQ_ini + 0.9074074074074074*hq0_raw*LQ_ini)