# Placeholder catalog

This directory stands in for the real PUGS halo catalog so the container builds
end-to-end before the 2048³ catalog is ready. It contains no Parquet files, so
`pugs.io.read_catalog` will raise `FileNotFoundError` against it.

Rebuild with a real catalog:

    ./build_container.sh -c /path/to/NUGS2048_catalog
