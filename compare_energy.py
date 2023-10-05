import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import sys
import subprocess
import re
import numpy as np
import findpairs
import MDAnalysis as mda
import argparse


def printUsage():
    print("TestGromacsForces.py requires at least one argument of the form folder=(folder)")
    print("folder should be the name (without extension) of the top/gro/mdp files to load")
    print("It also accepts additional arguments of the form argumentName=value:")
    print("platform:  A platform (Reference, CPU, OpenCL, or CUDA); default is platform=Reference")
    print("precision: The precision to use (single, mixed, or double); default is single")
    print()
    print("Example:   python TestGromacsForces.py folder=dhfr.pbc platform=CUDA")

# Set default values
folder = None
gmxPath = 'PATH/TO/GMX'
gromacsTopDir = './'
platformName = 'Reference'
precision = 'double'

# Parse the argument list
for arg in sys.argv[1:]:
    argSplit = arg.split("=")
    argName = argSplit[0]
    argValue = argSplit[1]
    if argName=='folder':
        folder = argValue
    elif argName=='gmxPath':
        gmxPath = argValue
    elif argName=='gromacsTopDir':
        gromacsTopDir = argValue
    elif argName=='platform':
        platformName = argValue
    elif argName=='precision':
        precision = argValue


# Ensure that the argument folder is provided
if folder is None:
    print("Error: Argument list must contain a folder")
    printUsage()
    sys.exit(1) 
print("folder =", folder)
print('gmxPath = ', gmxPath)
print('gromacsTopDir =', gromacsTopDir)
print("platformName =", platformName)
print("precision =", precision)
if platformName == 'OpenCL':
    properties = {'OpenCLPrecision':precision}
elif platformName == 'CUDA':
    properties = {'CudaPrecision':precision}
else:
    properties = dict()


# Compute forces with Gromacs
top_file_gromacs = folder+'/for_gromacs.top'
gro_file = folder+'/coords.gro'

subprocess.call([gmxPath, 'grompp', '-f', 'gmx_params.mdp', '-p', top_file_gromacs, '-c', gro_file, '-o', folder+'.tpr', '-maxwarn', '1'])
subprocess.call([gmxPath, 'mdrun', '-deffnm', folder, '-nt', '1'])
process = subprocess.Popen([gmxPath, 'dump', '-f', folder+'.trr'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(output, errors) = process.communicate()
expr = re.compile('f\[.*?\]=\{(.*?),(.*?),(.*?)\}')
forces_gromacs = []
for match in re.findall(expr, output.decode()):
    forces_gromacs.append(mm.Vec3(*[float(x) for x in match]))
terms = {}
terms['bonds'] = {'Bond', 'U-B'}
terms['angles'] = {'Angle'}
terms['torsions'] = {'Proper Dih.', 'Periodic Improper Dih.', 'CMAP Dih.'}
terms['nonbonded'] = {'LJ-14', 'Coulomb-14', 'LJ (SR)', 'Coulomb (SR)', 'Coul. recip.', 'Disper. corr.'}
energies_gromacs = dict((term, 0.0) for term in terms)
process = subprocess.Popen([gmxPath, 'dump', '-e', folder+'.edr'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(output, errors) = process.communicate()
for line in output.decode().split('\n'):
    line = line.strip()
    for term in terms:
        for prefix in terms[term]:
            if line.startswith(prefix):
                energies_gromacs[term] += float(line[len(prefix):].split()[0])


# Compute forces with OpenMM
top_file_gromacs= folder+'/for_gromacs.top'
gro_file = folder+'/coords.gro'
gro = app.GromacsGroFile(gro_file)
top = app.GromacsTopFile(top_file_gromacs, periodicBoxVectors=gro.getPeriodicBoxVectors(), includeDir=gromacsTopDir)
mdpOptions = {}
for line in open('gmx_params.mdp'):
    if '=' in line:
        fields = [x.strip().lower() for x in line.split('=')]
    mdpOptions[fields[0]] = fields[1]
if mdpOptions.get('implicit_solvent') == 'gbsa':
    system = top.createSystem(nonbondedMethod=app.NoCutoff, implicitSolvent=app.OBC2)
elif mdpOptions.get('coulombtype') == 'pme':
    system = top.createSystem(nonbondedMethod=app.PME, ewaldErrorTolerance=5e-5, nonbondedCutoff=0.9*unit.nanometers)
else:
    system = top.createSystem(nonbondedMethod=app.NoCutoff)
integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
simulation = app.Simulation(top.topology, system, integrator, mm.Platform.getPlatformByName(platformName), properties)
simulation.context.setPositions(gro.positions)
simulation.context.applyConstraints(1e-6)
state = simulation.context.getState(getForces=True)
forces_openmm = state.getForces().value_in_unit(unit.kilojoules_per_mole/unit.nanometer)
forces_openmm =np.array(forces_openmm)
terms = {}
terms['bonds'] = {'HarmonicBondForce'}
terms['angles'] = {'HarmonicAngleForce'}
terms['torsions'] = {'PeriodicTorsionForce'}
terms['nonbonded'] = {'NonbondedForce', 'CustomNonbondedForce'}
energies_openmm = dict((term, 0.0) for term in terms)
for i,f in enumerate(system.getForces()):
    f.setForceGroup(i)
for i,f in enumerate(system.getForces()):
    energy_i = simulation.context.getState(getEnergy=True, groups={i}).getPotentialEnergy()
    for term in terms:
        if f.getName() in terms[term]:
            energies_openmm[term] += energy_i._value


# Compute forces with modified OpenMM
top_file_modified_openmm= folder+'/for_modified_openmm.top'
gro_file = folder+'/coords.gro'
gro = app.GromacsGroFile(gro_file)
top = app.GromacsTopFile(top_file_modified_openmm, periodicBoxVectors=gro.getPeriodicBoxVectors(), includeDir=gromacsTopDir)
mdpOptions = {}
for line in open('gmx_params.mdp'):
    if '=' in line:
        fields = [x.strip().lower() for x in line.split('=')]
    mdpOptions[fields[0]] = fields[1]
if mdpOptions.get('implicit_solvent') == 'gbsa':
    system = top.createSystem(nonbondedMethod=app.NoCutoff, implicitSolvent=app.OBC2)
elif mdpOptions.get('coulombtype') == 'pme':
    system = top.createSystem(nonbondedMethod=app.PME, ewaldErrorTolerance=5e-5, nonbondedCutoff=0.9*unit.nanometers)
else:
    system = top.createSystem(nonbondedMethod=app.NoCutoff)
#with open('system_bf_modified.xml', 'w') as output:
#    output.write(mm.XmlSerializer.serialize(system))
## load topology and calculate special pairs that need to be modified
u = mda.Universe(top_file_modified_openmm, topology_format='ITP',infer_system=True, include_dir=gromacsTopDir)
atp1, atp2 = 'OB', 'HB'
combined_sigma = 0.150
combined_epsilon = 1.2552
pair_1_2 = findpairs.find_1_2_pairs(u, atp1, atp2)
pair_1_3 = findpairs.find_1_3_pairs(u, atp1, atp2)
pair_1_4 = findpairs.find_1_4_pairs(u, atp1, atp2)
pair_all = findpairs.find_all_pairs(u, atp1, atp2)
pair_SR = sorted(list(pair_1_2 | pair_1_3))
pair_LR = sorted(list(pair_all - pair_1_2 - pair_1_3 - pair_1_4))
pair_1_4 = sorted(list(pair_1_4))
param_SR = findpairs.get_nb_param(pair_SR, combined_sigma, combined_epsilon, u.atoms.charges, 'SR')
param_LR = findpairs.get_nb_param(pair_LR, combined_sigma, combined_epsilon, u.atoms.charges, 'LR')
param_1_4 = findpairs.get_nb_param(pair_1_4, combined_sigma, combined_epsilon, u.atoms.charges, '1-4')

for i,f in enumerate(system.getForces()):
    if isinstance(f, mm.NonbondedForce):
        nb_force = f   # find nb force

for pairs, params in [[pair_SR, param_SR], [pair_LR, param_LR], [pair_1_4, param_1_4]]:
    assert len(pairs) == len(params)
    for i, (at1, at2) in enumerate(pairs):
        eps_i, qq_i, sig_i = params[i]
        nb_force.addException(at1, at2, qq_i, sig_i, eps_i, replace=True)

#with open('system_af_modified.xml', 'w') as output:
#    output.write(mm.XmlSerializer.serialize(system))

integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
simulation = app.Simulation(top.topology, system, integrator, mm.Platform.getPlatformByName(platformName), properties)
simulation.context.setPositions(gro.positions)
simulation.context.applyConstraints(1e-6)
state = simulation.context.getState(getForces=True)
forces_openmm_modified = state.getForces().value_in_unit(unit.kilojoules_per_mole/unit.nanometer)
energies_openmm_modified = dict((term, 0.0) for term in terms)
for i,f in enumerate(system.getForces()):
    f.setForceGroup(i)
for i,f in enumerate(system.getForces()):
    energy_i = simulation.context.getState(getEnergy=True, groups={i}).getPotentialEnergy()
    for term in terms:
        if f.getName() in terms[term]:
            energies_openmm_modified[term] += energy_i._value

# Report results
print('###################### gromacs VS original openmm ######################')
print()

relativeDiff_1 = [2*np.linalg.norm(f1-f2)/(np.linalg.norm(f1)+np.linalg.norm(f2)) for f1, f2 in zip(forces_gromacs, forces_openmm)]
relativeDiff_1.sort()
print('Relative force error')
print('max:', relativeDiff_1[-1])
print('mean: ', sum(relativeDiff_1)/len(relativeDiff_1))
print('median:', relativeDiff_1[len(relativeDiff_1)//2])
print()

projection = [np.dot(f1, f2)/np.dot(f2, f2) for f1, f2 in zip(forces_gromacs, forces_openmm)]
projection = np.array(projection).reshape(-1,1)
print('Projection')
print('mean:', sum(projection)/len(projection))
print('min:', min(projection))
print('max:', max(projection))
print()

print('Potential energy')
print('\tgromacs (kj/mol)\topenmm (kj/mol)\trelative error')
for term in terms:
    e_gmx = energies_gromacs[term]
    e_openmm = energies_openmm[term]
    rela_diff = abs(e_gmx - e_openmm) / abs(e_gmx) if e_gmx != 0 else 0.0
    print(f'{e_gmx}\t{e_openmm}\t{rela_diff}')
print()
print('###################### gromacs VS modified openmm ######################')
print()

relativeDiff_1 = [2*np.linalg.norm(f1-f2)/(np.linalg.norm(f1)+np.linalg.norm(f2)) for f1, f2 in zip(forces_gromacs, forces_openmm_modified)]
relativeDiff_1.sort()
print('Relative force error')
print('max:', relativeDiff_1[-1])
print('mean: ', sum(relativeDiff_1)/len(relativeDiff_1))
print('median:', relativeDiff_1[len(relativeDiff_1)//2])
print()

projection = [np.dot(f1, f2)/np.dot(f2, f2) for f1, f2 in zip(forces_gromacs, forces_openmm_modified)]
projection = np.array(projection).reshape(-1,1)
print('Projection')
print('mean:', sum(projection)/len(projection))
print('min:', min(projection))
print('max:', max(projection))
print()

print('Potential energy')
print('\tgromacs (kj/mol)\topenmm (kj/mol)\trelative error')
for term in terms:
    e_gmx = energies_gromacs[term]
    e_openmm = energies_openmm_modified[term]
    rela_diff = abs(e_gmx - e_openmm) / abs(e_gmx) if e_gmx != 0 else 0.0
    print(f'{e_gmx}\t{e_openmm}\t{rela_diff}')
