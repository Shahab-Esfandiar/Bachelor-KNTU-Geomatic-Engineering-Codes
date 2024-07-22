SELECT SUM(ST_AREA(ST_INTERSECTION(CITY_REGIONS.GEOM, SOIL.GEOM)))
FROM CITY_REGIONS JOIN SOIL ON ST_INTERSECTS(CITY_REGIONS.GEOM, SOIL.GEOM)
WHERE SOIL.SOIL='IVS'