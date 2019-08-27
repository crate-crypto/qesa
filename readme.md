Base:IPA_no_zk
Second layer: ipa_almost_zk
Third layer: Qesa_inner


So IPA_no_zk is a base layer, ipa_almost_zk uses ipa_no_zk and qesa_inner uses ipa_almost_zk. Almost_zk refers to SZK
----------

Qesa_copy uses qesa_inner
Qesa_zk uses qesa_inner

--------
Base: lampa_no_zk
second layer: lpma_simple_zk


What are the uses of each of these components?


IPA_NO_ZK : Proves that the inner product of two vectors (`l` and `r`) is equal to `t` without zero knowledge with O(log n ) communications. 

IPA_ALMOST_ZK: Proving the same statement as IPA_NO_ZK. uses IPA_NO_ZK as a foundation but makes the protocol statistical zero knowledge

QESA_INNER: is a part of natural separation of QESA_ZK. uses IPA_ALMOST_ZK once 