import dadi
import numpy as np

def model_func(params, ns, pts):
	nu_1, nu_2, t1, nu11, nu12, m1_12, m1_21, nu12_1, nu12_2, t2, nu21, nu22, nu23, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	nu1_func = lambda t: nu_1 + (nu11 - nu_1) * (t / t1)
	nu2_func = lambda t: nu_2 + (nu12 - nu_2) * (t / t1)
	phi = dadi.Integration.two_pops(phi, xx, T=t1, nu1=nu1_func, nu2=nu2_func, m12=m1_12, m21=m1_21)
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	nu2_func = lambda t: nu12_1 + (nu22 - nu12_1) * (t / t2)
	nu3_func = lambda t: nu12_2 + (nu23 - nu12_2) * (t / t2)
	phi = dadi.Integration.three_pops(phi, xx, T=t2, nu1=nu21, nu2=nu2_func, nu3=nu3_func, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*len(ns))
	return sfs

dd = dadi.Misc.make_data_dict_vcf(vcf_filename='/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_STRICT_vcf_final.recode.vcf', popinfo_filename='/scratch/a_monc/postdoc/xipho_project/vcftools_filtering/xiph_pops.txt', filter=True, flanking_info=[None, None])
data = dadi.Spectrum.from_data_dict(dd, pop_ids=['Tap', 'Xin', 'Bel'], projections=[18, 24, 18], polarized=False)
pts = [30, 40, 50]
ns = data.sample_sizes

p0 = [0.010171591237009973, 0.01, 0.07626638644564733, 8.0141940197728, 1.6647910456390265, 9.17814477640171e-13, 6.4496383611851046e-12, 12.665498811125717, 0.010054534365352804, 0.024723805174405065, 16.237378589547223, 0.029646778614780316, 0.6300718758837248, 0.0, 0.0, 1.3584165336613545, 3.724689636449235, 0.6021256085909654, 6.279087993415007]
lower_bound = [0.01, 0.01, 1e-15, 0.01, 0.01, 0.0, 0.0, 0.01, 0.01, 1e-15, 0.01, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
upper_bound = [100.0, 100.0, 5.0, 100.0, 100.0, 10.0, 10.0, 100.0, 100.0, 5.0, 100.0, 100.0, 100.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
func_ex = dadi.Numerics.make_extrap_log_func(model_func)
model = func_ex(p0, ns, pts)
# we use fixed theta as we want to satisfy model limitations (e.g. bounds on time splits)
theta = 14966.98045600971
ll_model = dadi.Inference.ll(theta * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))

#Uncomment the next line to obtain the true optimal value of theta
#theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))

Nanc = 381150.8759664242
mu = 4.6e-09
L = 2134123
theta0 = 4 * mu * L
Nanc = int(theta / theta0)
print('Size of ancestral population: {0}'.format(Nanc))
