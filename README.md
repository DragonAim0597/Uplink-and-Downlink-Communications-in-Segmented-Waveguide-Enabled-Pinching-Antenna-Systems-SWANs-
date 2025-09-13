# Uplink and Downlink Communications in Segmented Waveguide-Enabled Pinching-Antenna Systems (SWANs)

The code for the paper **Z. Wang, C. Ouyang, and Y. Liu, “Uplink and downlink communications in segmented waveguide-enabled pinching-antenna systems (SWANs),” Sep. 2025.** 

<img decoding="async" src="./img/system_model.jpg" width="50%">

Abstract: A segmented waveguide-enabled pinching-antenna system (SWAN) is proposed, in which a segmented waveguide composed of multiple short dielectric waveguide segments is employed to radiate or receive signals through the pinching antennas (PAs) deployed on each segment. Based on this architecture, three practical operating protocols are proposed: segment selection (SS), segment aggregation (SA), and segment multiplexing (SM). For uplink SWAN communications, where one PA is activated per segment, the segmented structure eliminates the inter-antenna radiation effect, i.e., signals captured by one PA may re-radiate through other PAs along the same waveguide. This yields a tractable and physically consistent uplink signal model for a multi-PA pinching-antenna system (PASS), which has not been established for conventional PASS using a single long waveguide. Building on this model, PA placement algorithms are proposed to maximize the uplink signal-to-noise ratio (SNR). Closed-form expressions for the received SNR under the three protocols are derived, and the corresponding scaling laws with respect to the number of segments are analyzed. It is proven that the segmented architecture reduces both the average PA-to-user distance and the PA-to-feed distance, thereby mitigating both large-scale path loss and in-waveguide propagation loss. These results are extended to downlink SWAN communications, where multiple PAs are activated per segment, and PA placement methods are proposed to maximize the downlink received SNR under the three protocols. Numerical results demonstrate that: i) among the three protocols, SM achieves the best performance, followed by SA and then SS; and ii) for all protocols, the proposed SWAN achieves a higher SNR than conventional PASS with a single long waveguide in both uplink and downlink scenarios.


## Running the simulations

### Prerequisites

- [MATLAB](https://uk.mathworks.com/products/matlab.html)

### Launch

Run `Uplink_SWAN.m` to plot Fig. 9 in this paper.

Run `Downlink_SWAN.m` to plot Fig. 11 in this paper.

```


