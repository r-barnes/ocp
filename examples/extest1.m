clear all
pc6 = load_ocp("extest1const6.data"); 
pc11 = load_ocp("extest1const11.data");
pc21 = load_ocp("extest1const21.data");
pc41 = load_ocp("extest1const41.data");
pc81 = load_ocp("extest1const81.data");

printf("CONSTANT\n");
extest1_max_err(pc6)
extest1_max_err(pc11)
extest1_max_err(pc21)
extest1_max_err(pc41)
extest1_max_err(pc81)

pl6 = load_ocp("extest1lin6.data"); 
pl11 = load_ocp("extest1lin11.data");
pl21 = load_ocp("extest1lin21.data");
pl41 = load_ocp("extest1lin41.data");
pl81 = load_ocp("extest1lin81.data");

printf("LINEAR\n");
extest1_max_err(pl6)
extest1_max_err(pl11)
extest1_max_err(pl21)
extest1_max_err(pl41)
extest1_max_err(pl81)

pu6 = load_ocp("extest1cubic6.data"); 
pu11 = load_ocp("extest1cubic11.data");
pu21 = load_ocp("extest1cubic21.data");

printf("CUBIC\n");
extest1_max_err(pu6)
extest1_max_err(pu11)
extest1_max_err(pu21)
