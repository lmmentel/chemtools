# coding: utf-8
import chemtools.gamessus as gus
gg = gus.Gamess(
        name='GamessUS',
        execpath='/home/lmentel/Programs/gamess-us-may2013',
        scratch='/home/lmentel/scratch',
        )
print "Attributes:"
print "\t{0:<15s} : {1:}".format("name", gg.name)
print "\t{0:<15s} : {1:}".format("execpath", gg.execpath)
print "\t{0:<15s} : {1:}".format("runopts", gg.runopts)
print "\t{0:<15s} : {1:}".format("scratch", gg.scratch)
print "\t{0:<15s} : {1:}".format("rungms", gg.rungms)
print "\t{0:<15s} : {1:}".format("version", gg.version)
print "\t{0:<15s} : {1:}".format("ddikick", gg.ddikick)

gg.run(inpfile="ne_dz_guga.inp")

log = gus.GamessLogParser("ne_dz_guga.log")

print "{0:<30s} : {1:}".format("terminatedOK", log.accomplished())
print "{0:<30s} : {1:}".format("gamess version", log.get_version())
print "{0:<30s} : {1:}".format("charge", log.get_charge())
print "{0:<30s} : {1:}".format("electrons", log.get_electrons())
print "{0:<30s} : {1:}".format("homo", log.get_homo())
print "{0:<30s} : {1:}".format("number of atoms", log.get_number_of_atoms())
print "{0:<30s} : {1:}".format("number of AOs", log.get_number_of_aos())
print "{0:<30s} : {1:}".format("number of MOs", log.get_number_of_mos())
print "{0:<30s} : {1:}".format("linear deps", log.get_linear_deps())
print "{0:<30s} : {1:}".format("SCF type", log.get_scf_type())
print "{0:<30s} : {1:}".format("CC type", log.get_cc_type())
print "{0:<30s} : {1:}".format("CI type", log.get_ci_type())

print "\n Components of the HF energy"
for key, value in log.get_energy_components("hf").items():
    print "{0:<40s} = {1:>18.10f}".format(key, value)

print "\n Components of the CI energy"
for key, value in log.get_energy_components("ci").items():
    print "{0:<40s} = {1:>18.10f}".format(key, value)
