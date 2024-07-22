SELECT SUM(BALANCE_ACC)
FROM ((public.depositer JOIN public.account ON ACC_DEP=ACC_ID) JOIN public.customer ON NAME_DEP=NAME_CUS)
WHERE CITY_CUS='Isfahan'