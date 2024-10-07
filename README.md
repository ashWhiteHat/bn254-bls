# BLS Signature on Bn254
BLS signature implementation on bn254 (alt_bn128 / bn256) curve according to [EIP-196](https://eips.ethereum.org/EIPS/eip-196).

$G1: y^2 = x^3 + 3$

$G2: y^2 = x^3 + 3(u + 1)$

These two group supports bilinearity by pairing. Let $G$ and $H$ be generator of $G1$, and $G2$, and $e$ be pairing function. The relationship is described as following.

$e(aG, bH) = e(G, H)^{ab}$

## Test

```shell
$ cargo test
```
