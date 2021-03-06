using Images

function init_O(I, H_high, W_high, mode="wfp")
    H_low, W_low, L = size(I)
        
    if mode == "wfp"
        o_im = sqrt.( Complex.(I[:,:,trunc(Int,L/2)+1] )) * (H_low*W_low/H_high/W_high)
        o_im = imresize(o_im, H_high, W_high)
        O0 = fftshift(fft(o_im))
    elseif mode == "gs"
        A = imresize(I[:,:,trunc(Int,L/2)+1], H_high, W_high)
        Mag_image = H_high / H_low
        scale = Mag_image.^2
        A = A./scale
        Hi_res_M, Hi_res_N = size(A)
        F = fftshift(fft(A))
    end
end