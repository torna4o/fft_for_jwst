# Convenience for several FFT filters for JWST MIRI imaging

#####################################
#####################################
#                                   #
#                                   #
#    Loading Required Libraries     #
#                                   #
#                                   #
#####################################
#####################################

import numpy as np # Numerical Python, fft and other numpy array operations
import matplotlib.pyplot as plt # Plotting figures

#####################################
#####################################
#                                   #
#                                   #
#            Functions              #
#                                   #
#                                   #
#####################################
#####################################

def fft_circmask_pad(dat, threshold, plotting=True, vm=0, vm2=0, vma=1, vma2=1):
    
    """
    This function is special for padded input data
    The input is not square in 2D, so it previously
    had padding with zero till reaching a two-power
    dimensions.
    Then, after FFT, it returns back to its before-pad
    shape.
    """
    
    
    """
    dat is data object, a numpy.ndarray, in this case
        zero-padded to 1024x1024 pixels
        
    threshold is the radius of the circle that will 
        filter out a part of data within the spectral domain
        the bigger, the more data will be lost
        
    plotting is a switch, put True if you would like to
        plot the FFT filtered and the filtrate images
        side-by-side after the operation
    
    vm, vm2 and vma, vma2 are minimum and maximum values to be
        used in matplotlib.pyplot.imshow() plotting part.
        vm and vma are for original and filtered images, while
        vm2 and vma2 are for filtrate.
        The point is, if you have a filter that leaves low amount
        of data in the filtrate, you would like to have low values
        for vm2 and vma2.
    """
    img_fft = np.fft.fft2(dat)
    img_fft2 = np.fft.fftshift(img_fft) # Shift the zero-freq. to center
    
    mask = np.zeros(dat.shape)
    mask_inv = np.ones(dat.shape)
    center = (mask.shape[0] // 2, mask.shape[1] // 2)
    
    """
    The following np.ogrid creates an index grid
    that will allow pointing to specific
    locations later
    """

    y, x = np.ogrid[:mask.shape[0], :mask.shape[1]]
    
    """
    Now we obtain mask after array operation
    More precisely, "mask_area" is a True False matrix, 
    it will show whether a specific indexed data will
    be included in the final masked numpy array.
    """

    mask_area = (x - center[1])**2 + (y - center[0])**2 <= threshold**2
    mask[mask_area] = 1
    mask_inv[mask_area] = 0
    
    # Apply the mask to the FFT output
    img_fft_filtered = img_fft2 * mask
    img_fft_filtered_inv = img_fft2 * mask_inv
    
    # Shift the zero-frequency component back to the top-left corner of the spectrum
    img_fft_filtered_k = np.fft.ifftshift(img_fft_filtered)
    img_fft_filtered_inv_k = np.fft.ifftshift(img_fft_filtered_inv)

    # Apply inverse FFT to get the filtered image
    img_filtered = np.abs(np.fft.ifft2(img_fft_filtered_k)) # np.abs calculates absolute values
    img_filtered_inv = np.abs(np.fft.ifft2(img_fft_filtered_inv_k)) # np.abs calculates absolute values

    img_filtered_z = img_filtered[50:(50 + 780), 50:(50 + 560)]
    img_filtered_inv = img_filtered_inv[50:(50 + 780), 50:(50 + 560)]
    print(img_filtered.shape)
    print("### *** \\\ FFT TASK COMPLETE //// **** ###")
    
    if plotting==True:
       
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(2, 3)
        axs[0, 0].imshow(dat[50:(50 + 780), 50:(50 + 560)], origin='lower', cmap='afmhot', vmin= vm, vmax=vma)
        axs[0, 0].set_title('Image before FFT')
        axs[0, 1].imshow(img_filtered_z, origin='lower', cmap='afmhot', vmin= vm, vmax= vma)
        axs[0, 1].set_title('Filtered Image')
        axs[0, 2].imshow(img_filtered_inv, origin='lower', cmap='afmhot', vmin= vm2, vmax= vma2)
        axs[0, 2].set_title('The Filtrate')
        axs[1, 0].imshow((np.abs(img_fft2)), cmap='gray', vmin= 0, vmax=8000)
        axs[1, 0].set_title('Original Image FFT')
        axs[1, 1].imshow(np.abs(img_fft_filtered), cmap='gray', vmin=0, vmax= 8000)
        axs[1, 1].set_title('Mask FFT')
        axs[1, 2].imshow(np.abs(img_fft_filtered_inv), cmap='gray', vmin=0, vmax= 8000)
        axs[1, 2].set_title('Inverse of the Mask for FFT')
        plt.show()
            
    return img_filtered, img_filtered_inv


def fft_linemask_pad(dat, vm, vm2, vma1, vma2, v_thick = 0, h_thick = 0, plotting=True):

    """
    This function is special for padded input data
    The input is not square in 2D, so it previously
    had padding with zero till reaching a two-power
    dimensions.
    Then, after FFT, it returns back to its before-pad
    shape.
    """
    
    
    """
    dat is data object, a numpy.ndarray, in this case
        zero-padded to 1024x1024 pixels
        
    v_thick is thickness of the vertical line to be removed
        starting from the center of the spectral dimension
        it can be set 0 to not removing anything in that axis
    
    h_thick is the horizontal counterpart of the v_thick
    
    plotting is a switch, put True if you would like to
        plot the FFT filtered and the filtrate images
        side-by-side after the operation
    vm, vm2 and vma, vma2 are minimum and maximum values to be
        used in matplotlib.pyplot.imshow() plotting part.
        vm and vma are for original and filtered images, while
        vm2 and vma2 are for filtrate.
        The point is, if you have a filter that leaves low amount
        of data in the filtrate, you would like to have low values
        for vm2 and vma2.
    """
    img_fft = np.fft.fft2(dat)
    img_fft_filtered = np.fft.fftshift(img_fft) # Shift the zero-freq. to center
    img_fft_filtered_back = np.fft.fftshift(img_fft) # Shift the zero-freq. to center

    img_fft_filtered_inv_sh = np.zeros((img_fft.shape[0], img_fft.shape[1]))
    print("Shape of Zero matrix for inv:", img_fft_filtered_inv_sh.shape)
        
    if v_thick > 0:
        img_fft_filtered_inv_sh[:(img_fft.shape[0])//2, (img_fft.shape[1])//2-v_thick:(img_fft.shape[1])//2+v_thick+1] = np.abs(img_fft[:(img_fft.shape[0])//2, (img_fft.shape[1])//2-v_thick:(img_fft.shape[1])//2+v_thick+1])
        img_fft_filtered_inv_sh[-img_fft.shape[0]//2:, (img_fft.shape[1]//2)-v_thick:(img_fft.shape[1]//2)+v_thick+1] = np.abs(img_fft[-img_fft.shape[0]//2:, (img_fft.shape[1]//2)-v_thick:(img_fft.shape[1]//2)+v_thick+1])
        img_fft_filtered[:img_fft.shape[0]//2, (img_fft.shape[1]//2)-v_thick:(img_fft.shape[1]//2)+v_thick+1] = 0
        img_fft_filtered[-img_fft.shape[0]//2:, (img_fft.shape[1]//2)-v_thick:(img_fft.shape[1]//2)+v_thick+1] = 0
    if h_thick > 0:
        img_fft_filtered_inv_sh[(img_fft.shape[0]//2)-h_thick:(img_fft.shape[0]//2)+h_thick+1,:img_fft.shape[1]//2] = np.abs(img_fft[(img_fft.shape[0]//2)-h_thick:(img_fft.shape[0]//2)+h_thick+1,:img_fft.shape[1]//2])
        img_fft_filtered_inv_sh[(img_fft.shape[0]//2)-h_thick:(img_fft.shape[0]//2)+h_thick+1,-img_fft.shape[1]//2:] = np.abs(img_fft[(img_fft.shape[0]//2)-h_thick:(img_fft.shape[0]//2)+h_thick+1,-img_fft.shape[1]//2:])
        img_fft_filtered[(img_fft.shape[0]//2)-h_thick:(img_fft.shape[0]//2)+h_thick+1,:img_fft.shape[1]//2] = 0
        img_fft_filtered[(img_fft.shape[0]//2)-h_thick:(img_fft.shape[0]//2)+h_thick+1,-img_fft.shape[1]//2:] = 0


    # Shift the zero-frequency component back to the top-left corner of the spectrum
    img_fft_filtered2 = np.fft.ifftshift(img_fft_filtered)
    img_fft_filtered_inv = np.fft.ifftshift(img_fft_filtered_inv_sh)

    # Apply inverse FFT to get the filtered image
    img_filtered = np.abs(np.fft.ifft2(img_fft_filtered2)) # np.abs calculates absolute values
    img_filtered_inv_inv = np.abs(np.fft.ifft2(img_fft_filtered_inv)) # np.abs calculates absolute values

    img_filtered_k = img_filtered[50:(50 + 780), 50:(50 + 560)]
    img_filtered_inv_inv = img_filtered_inv_inv[50:(50 + 780), 50:(50 + 560)]
    print(img_filtered.shape)
    print("### *** \\\ FFT TASK COMPLETE //// **** ###")
    
    if plotting==True:
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(2, 3)
        axs[0, 0].imshow(dat[50:(50 + 780), 50:(50 + 560)], origin='lower', cmap='afmhot', vmin= vm, vmax=vma)
        axs[0, 0].set_title('Image before FFT')         
        axs[0, 1].imshow(img_filtered_k, origin='lower', cmap='afmhot', vmin= vm, vmax= vma)
        axs[0, 1].set_title('Filtered Image')
        axs[0, 2].imshow(img_filtered_inv_inv, origin='lower', cmap='afmhot', vmin= vm2, vmax= vma2)
        axs[0, 2].set_title('The Filtrate')
        axs[1, 0].imshow(np.abs(img_fft_filtered_back), cmap='gray', vmin= 0, vmax=200)
        axs[1, 0].set_title('Original Image FFT')
        axs[1, 1].imshow(np.abs(img_fft_filtered), cmap='gray', vmin=4, vmax=8000)
        axs[1, 1].set_title('Mask FFT')
        axs[1, 2].imshow(np.abs(img_fft_filtered_inv_sh), cmap='gray', vmin=0, vmax= 100)
        axs[1, 2].set_title('Inverse Mask FFT')        
        
        plt.show()        

    return img_filtered, img_filtered_inv_inv

#####################################
#####################################
#                                   #
#                                   #
#            Execution              #
#                                   #
#                                   #
#####################################
#####################################


#aaa, bbb = fft_circmask_pad(dat=zero_padded, threshold=420, plotting=True, fitsingest=False, vm=21, vma=26, vm2=0.01, vma2=0.2)
#aaa, bbb = fft_linemask_pad(dat=zero_padded, v_thick=10, h_thick=10, plotting=True, p_vmin1=5.95, p_vmin2=0.06, p_vmax1=6.8, p_vmax2=4)
