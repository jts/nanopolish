# Developer Notes

## Updating Bioconda on tagged releases
The following is a quick step-by-step checklist on updating the bioconda release for nanopolish, to be done after each tagged release, and is a condensed/updated version of [these slides](https://monashbioinformaticsplatform.github.io/bioconda-tutorial/#/) by Andrew Perry.
1. On Github, fork `https://github.com/bioconda/bioconda-recipes` to `https://github.com/{USER}/bioconda-recipes` and clone the latter repository to a local directory; `cd` into the cloned directory.
2. Check out a new branch via `git branch nanopolish-bioconda-bump && git checkout nanopolish-bioconda-bump`.
3. Update the `bioconda-recipes/recipes/nanopolish/meta.yaml` file by editing the version tag and the SHA hash; the SHA256 hash can be obtained by running `sha256sum nanopolish-v{$VERSION}.tar.gz` on the command line (where `{VERSION}` is the new, updated version tag); commit the changes to the `meta.yaml` file via, e.g., `git commit -a -m 'bump nanopolish to version {VERSION}'`.
4. Push the changes to your forked repo via `git push origin nanopolish-bioconda-bump`; then, make a pull request to merge the updates into the master branch of the upstream `bioconda-recipes` repository.
5. If all goes well, the automated TravisCI tests on the upstream repository will pass and an owner will merge the changes.
6. Otherwise, if further edits are requested or if the TravisCI tests fail, make further commits to the local cloned repository and push to the forked repository on Github; the changes should automatically appear in the pull request and will trigger an automated TravisCI check.
