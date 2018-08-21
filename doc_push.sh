#!/bin/bash
################################################################################
# Modified from https://gist.github.com/vidavidorra/548ffbcdae99d752da02

################################################################################
##### Setup this script and get the current gh-pages branch.               #####
echo 'Setting up the script...'
# Exit with nonzero exit code if anything fails
set -e

# get path to Doxyfile
REPOPATH=$(pwd)
DOXYFILE="$REPOPATH/Doxyfile"

# set other variables
GH_REPO_NAME=accelerInt
GH_REPO_REF=git@github.com:SLACKHA/accelerInt.git
GH_COMMIT_FULL="$(git rev-parse HEAD)"
#take substring to get short hash
GH_COMMIT=${GH_COMMIT_FULL:0:6}

# Create a clean working directory for this script.
mkdir -p code_docs
cd code_docs

# Get the current gh-pages branch
git clone -b gh-pages $GH_REPO_REF
cd $GH_REPO_NAME

##### Configure git.
# Set the push default to simple i.e. push only the current branch.
git config --global push.default simple

# Remove everything currently in the gh-pages branch.
# GitHub is smart enough to know which files have changed and which files have
# stayed the same and will only update the changed files. So the gh-pages branch
# can be safely cleaned, and it is sure that everything pushed later is the new
# documentation.
rm -rf *

# Need to create a .nojekyll file to allow filenames starting with an underscore
# to be seen on the gh-pages site. Therefore creating an empty .nojekyll file.
# Presumably this is only needed when the SHORT_NAMES option in Doxygen is set
# to NO, which it is by default. So creating the file just in case.
echo "" > .nojekyll

################################################################################
##### Generate the Doxygen code documentation and log the output.          #####
echo 'Generating Doxygen code documentation...'
# Redirect both stderr and stdout to the log file AND the console.
REPOPATH=$REPOPATH doxygen $DOXYFILE 2>&1 | tee doxygen.log

################################################################################
##### Upload the documentation to the gh-pages branch of the repository.   #####
# Only upload if Doxygen successfully created the documentation.
# Check this by verifying that the file index.html
# both exist. This is a good indication that Doxygen did it's work.
if [ -f "index.html" ]; then

    echo 'Uploading documentation to the gh-pages branch...'
    # Add everything in this directory (the Doxygen code documentation) to the
    # gh-pages branch.
    # GitHub is smart enough to know which files have changed and which files have
    # stayed the same and will only update the changed files.
    git add --all

    # Commit the added files with a title and description containing the Travis CI
    # build number and the GitHub commit reference that issued this build.
    git commit -m "Deploy code docs to GitHub Pages for commit ${GH_COMMIT}".

    # Force push to the remote gh-pages branch.
    # The ouput is redirected to /dev/null to hide any sensitive credential data
    # that might otherwise be exposed.
    git push --force "${GH_REPO_REF}" > /dev/null 2>&1
else
    echo '' >&2
    echo 'Warning: No documentation (html) files have been found!' >&2
    echo 'Warning: Not going to push the documentation to GitHub!' >&2
    exit 1
fi