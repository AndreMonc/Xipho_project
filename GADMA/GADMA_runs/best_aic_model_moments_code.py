import moments
import numpy as np

def model_func(params, ns):
	nu_1, nu_2, t1, nu11, nu12, m1_12, m1_21, nu12_1, nu12_2, t2, nu21, nu22, nu23, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	sts = moments.LinearSystem_1D.steady_state_1D(np.sum(ns))
	fs = moments.Spectrum(sts)
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
	nu1_func = lambda t: nu_1 + (nu11 - nu_1) * (t / t1)
	nu2_func = lambda t: nu_2 + (nu12 - nu_2) * (t / t1)
	migs = np.array([[0, m1_12], [m1_21, 0]])
	fs.integrate(tf=t1, Npop=lambda t: [nu1_func(t), nu2_func(t)], m=migs, dt_fac=0.01)
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	nu2_func = lambda t: nu12_1 + (nu22 - nu12_1) * (t / t2)
	nu3_func = lambda t: nu12_2 + (nu23 - nu12_2) * (t / t2)
	migs = np.array([[0, m2_12, m2_13], [m2_21, 0, m2_23], [m2_31, m2_32, 0]])
	fs.integrate(tf=t2, Npop=lambda t: [nu21, nu2_func(t), nu3_func(t)], m=migs, dt_fac=0.01)
	return fs

dd = moments.Misc.make_data_dict_vcf(vcf_filename='/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_STRICT_vcf_final.recode.vcf', popinfo_filename='/scratch/a_monc/postdoc/xipho_project/vcftools_filtering/xiph_pops.txt', filter=True, flanking_info=[None, None])
data = moments.Spectrum.from_data_dict(dd, pop_ids=['Tap', 'Xin', 'Bel'], projections=[18, 24, 18], polarized=False)
ns = data.sample_sizes

p0 = [0.010171591237009973, 0.01, 0.07626638644564733, 8.0141940197728, 1.6647910456390265, 9.17814477640171e-13, 6.4496383611851046e-12, 12.665498811125717, 0.010054534365352804, 0.024723805174405065, 16.237378589547223, 0.029646778614780316, 0.6300718758837248, 0.0, 0.0, 1.3584165336613545, 3.724689636449235, 0.6021256085909654, 6.279087993415007]
lower_bound = [0.01, 0.01, 1e-15, 0.01, 0.01, 0.0, 0.0, 0.01, 0.01, 1e-15, 0.01, 0.01, 0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
upper_bound = [100.0, 100.0, 5.0, 100.0, 100.0, 10.0, 10.0, 100.0, 100.0, 5.0, 100.0, 100.0, 100.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
model = model_func(p0, ns)
# we use fixed theta as we want to satisfy model limitations (e.g. bounds on time splits)
theta = 14966.98045600971
ll_model = moments.Inference.ll(theta * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))

#Uncomment the next line to obtain the true optimal value of theta
#theta = moments.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))

Nanc = 381150.8759664242
mu = 4.6e-09
L = 2134123
theta0 = 4 * mu * L
Nanc = int(theta / theta0)
print('Size of ancestral population: {0}'.format(Nanc))


plot_ns = [4 for _ in ns]  # small sizes for fast drawing
gen_mod = moments.ModelPlot.generate_model(model_func,
                                           p0, plot_ns)
moments.ModelPlot.plot_model(gen_mod,
                             save_file='model_from_GADMA.png',
                             fig_title='Demographic model from GADMA',
                             draw_scale=True,
                             pop_labels=['Tap', 'Xin', 'Bel'],
                             nref=381150,
                             gen_time=0.003,
                             gen_time_units='kya',
                             reverse_timeline=True)