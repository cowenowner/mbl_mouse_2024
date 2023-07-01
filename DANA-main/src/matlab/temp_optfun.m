function o = temp_optfun(ISIs, target_lv)
lv = LocalVariance(ISIs);
o = lv-target_lv;