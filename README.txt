===========
Papillon
===========

A Python module to read cummeRbund generated files, mainly after
RNA-seq analysis through Galaxy Platform https://usegalaxy.org/

Typical usage:

    #!/usr/bin/env python

    import papillon as pp

    YourExperiment=pp.read_db("Your Folder")
    YourExperiment.get_gene()
    YourExperiment.heatmap()
