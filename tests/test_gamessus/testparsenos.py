# coding: utf-8

from chemtools.gamessus import GamessDatParser
gdp = GamessDatParser(datfile="ne_dz_guga.dat")
nos = gdp.get_nos()
orbs = gdp.parse_orbitals(nos)
print orbs[:, 0]
print orbs[:, -1]
#gdp = GamessDatParser(datfile="Yb_CPC_DK3_ANO-RCC-MB84_CCSDT-48e_06.00.dat")
#mos = gdp.get_nos()
