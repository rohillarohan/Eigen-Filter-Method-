# Eigenfilter-Method-
Eigen Filter Method is a method used to design filter with an arbitrary frequency response.

This method involves the computation of filter coefficients as the eigenvector of an appropriate Hermitian matrix. As opposed to the Least Squares approach where we need to find the matrix inverse to model a filter with the given frequency response, eigenfilter method only requires the computation of a single eigenvector. It is found that this method of filter design has low complexity as compared to the other methods as well it easily incorporates various time and frequency domain constraints. This method as numerous applications such as beamforming and channel shortening equalizers for discrete multitome systems.

This project has been divided in to two parts. 

In the first part I used eigenfilter approach to design a type-1 linear phase filter with 50 taps and sampling frequency of 24000 Hz whose frequency response should match the desired frequency response D(w).

In the second part I am essentially designing the same filter but this time the response have been made a little more arbitrary by including two notches at 4500 Hz and 8000 Hz.

At the end both frequency responses are compared with the desired response.
