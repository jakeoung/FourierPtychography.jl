using Plots
using FileIO
using TestImages

using FourierPtychography

# include("make_synthetic.jl")
# include("operator.jl")

img_real = Float64.(TestImages.testimage("cameraman"))
phase_temp = Float64.(TestImages.testimage("fabio_gray_512"))
phase_temp = (phase_temp .- minimum(phase_temp))
phase_temp ./= maximum(phase_temp)
phase_real = phase_temp .* (pi/2.0)

ratio_LR = 0.1
# ratio_step = ratio_LR * 0.4
ratio_step = ratio_LR * 0.4
sigmaN = 0.0

I, pupil, kxky, img, gt = make_synthetic(img_real, phase_real, ratio_LR, ratio_step, sigmaN)

# algorithm params

H_high, W_high = size(gt, 1), size(gt, 2)
O0 = init_O(I, H_high, W_high, "gs")
dresult = "result_gs/"
niter = 20
O, img_real_est, relerrs = recon_gs(I, O0, kxky, pupil, niter, dresult, gt)
