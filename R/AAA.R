## <FIXME> This is needed only for R 1.9.x, remove eventually

if(R.version$major < "2"){
    .First.lib <- function(...)
    {
        require("graphics")
        require("methods")
        require("stats")
        require("stats4")
    }

    .First.lib()
}
## </FIXME>
