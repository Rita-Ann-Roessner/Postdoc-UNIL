#!/bin/bash

find . -name "*sample*" | head -30
# To check files that will be deleted before doing it.
find . -name "*sample*" -delete
find . \( \! -name "*summary*" \) -name "*.json" | head -30
# Again check which files will be deleted.
find . \( \! -name "*summary*" \) -name "*.json" -delete
find . -name "TERMS_OF_USE.md" | head -30
find . -name "TERMS_OF_USE.md" -delete