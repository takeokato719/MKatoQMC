import subprocess

delta = '0.1'
# alpha = '0.1'
nsite = 256
nthermal = 1000
nbin = 100
nmeasure = 10
nseed = 1221

alpha_list = [ '0.1', '0.2', '0.3', '0.5', '1.0']

for alpha in alpha_list:
    runqmc = './runqmc nthermal %i nbin %i nmeasure %i nsite %i delta %s alpha %s nseed %i' % ( nthermal, nbin, nmeasure, nsite, delta, alpha, nseed )
    filename = 'corr_nthermal%i_nbin%i_nmeasure%i_nsite%i_delta%s_alpha%s_nseed%i.dat' % ( nthermal, nbin, nmeasure, nsite, delta, alpha, nseed )
    cmd = '%s > %s' % (runqmc, filename)
    print cmd
    subprocess.call(cmd, shell=True)


