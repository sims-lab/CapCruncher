Glossary
########

.. glossary::

    fragment
        The entire fragment of DNA that is sequenced. e.g.::

            --R1-->
            ----------fragment-----
                           <----R2--

    slice
        The results of *in silico* restriction digest of a fragment. e.g.::

        ---------GATC--------GATC-------------
        --slice1-    -slice2-    -slice3------

    capture
        A slice that overlaps one of the supplied viewpoints. e.g.::

            ---capture-slice-----
            -------viewpoint-----GATC--------------

    reporter
        A mapped slice from the same fragment as a capture slice used to report an interaction between
        two the viewpoint and another genomic region. e.g.::

            ---capture-----GATC--reporter-----
