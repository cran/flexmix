
### check potential coding problems via `codetools'
if (require("codetools")) {
    library("flexmix")
    library("nnet")
    library("mvtnorm")
    library("grid")
    print(checkUsagePackage("flexmix"))
}
