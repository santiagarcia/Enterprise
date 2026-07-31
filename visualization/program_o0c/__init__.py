"""O0C energy visuals.

A package rather than a bare directory so the visual gate can import the
renderer modules by name and scan their abstract syntax trees. The gate's whole
enforcement of "Python computes no physics" depends on reaching these modules
the same way in every configuration.
"""
