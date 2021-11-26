const BASE_URL = "https://km3netdbweb.in2p3.fr"
const COOKIE_FILENAME = expanduser("~/.km3netdb_cookie")
const SESSION_COOKIES = Dict(
    :lyon => "_kmcprod_134.158_lyo7783844001343100343mcprod1223user",
    :jupyter => "_jupyter-km3net_131.188.161.143_d9fe89a1568a49a5ac03bdf15d93d799",
    :gitlab => "_gitlab-km3net_131.188.161.155_f835d56ca6d946efb38324d59e040761",
)

const COOKIE_NAME = "sid"
const COOKIE_JAR = Dict{String, String}()

# _cookie_sid_pattern = re.compile(r"_[a-z0-9-]+_(\d{1,3}.){1,3}\d{1,3}_[a-z0-9]+")


function streamdsstreams()
    ascii2dataframe(request("streamds"))
end

function streamds(stream; kwargs...)
    ascii_data = request("streamds/$(stream).txt"; kwargs...)
    ascii2dataframe(ascii_data)
end

function detx(detid, run)
    runinfo = streamds("runs"; detid=detid, run=run)
    nrow(runinfo) == 0 && error("Run $(run) not found for det ID $(detid)")
    tcal = runinfo[1, :T0_CALIBSETID]
    pcal = runinfo[1, :POS_CALIBSETID]
    rcal = runinfo[1, :ROT_CALIBSETID]
    if tcal == ""
        @warn "Missing time calibration for det ID $(detid) run $(run)"
        tcal = 0
    end
    if pcal == ""
        @warn "Missing position calibration for det ID $(detid) run $(run)"
        pcal = 0
    end
    if rcal == ""
        @warn "Missing rotation calibration for det ID $(detid) run $(run)"
        rcal = 0
    end
    request("detx/$(detid)"; tcal=tcal, pcal=pcal, rcal=rcal)
end

function ascii2dataframe(ascii_data)
    !istabular(ascii_data) && error(ascii_data)
    CSV.read(IOBuffer(ascii_data), DataFrame)
end


"""
    function istabular(ascii_data; separator='\t')

Sanity check to verify that the ASCII data is tabular since unfortunately the
database API does no error reporting. The only way to see if a response is
an error is to check for line-breaks and tabs...
"""
function istabular(ascii_data; separator='\t')
    n_chars = length(ascii_data)
    n_chars == 0 && return false

    chunk = ascii_data[1:min(n_chars, 5000)]  # just an arbitrary large value
    lines = split(chunk, "\n")

    length(lines) < 2 && return false  # we require at least a header (or two data rows)

    separators_in_header = count(c->c==separator, lines[1])

    # two lines with equal amount of tabs is enough to be sure that it's tabulated
    separators_in_header > 0 && separators_in_header == count(c->c==separator, lines[2]) && return true

    false
end


function request(url; kwargs...)
    parameters = ""
    if length(kwargs) > 0
        parameters = "?"
        for (key, value) ∈ kwargs
            parameters *= "$(key)=$(value)&"
        end
    end

    url = "$(BASE_URL)/$(url)$(parameters)"
    r = HTTP.request("GET", url; cookies=COOKIE_JAR)
    r.status == 200 || error(r.error)
    length(r.body) == 0 && error("No data returned for $(url)")
    String(r.body)
end

function initdb()
    cookie = acquire_session_cookie()
    COOKIE_JAR[COOKIE_NAME] = cookie
    write_session_cookie(cookie)
end


function write_session_cookie(cookie)
    open(COOKIE_FILENAME, "w") do fobj
        write(fobj, ".in2p3.fr\tTRUE\t/\tTRUE\t0\tsid\t$(cookie)")
    end
    nothing
end


function acquire_session_cookie()
    haskey(ENV, "KM3NET_DB_COOKIE") && return ENV["KM3NET_DB_COOKIE"]

    if isfile(COOKIE_FILENAME)
        cookie = open(COOKIE_FILENAME) do f
            split(String(read(f)))[end]
        end
        return cookie
    end

    return cookie_via_login()
end


function cookie_via_login()
    print("Enter your KM3NeT DB username: ")
    user = readline()
    secret = Base.getpass("Password")
    password = String(read(secret))
    Base.shred!(secret)

    r = HTTP.request("GET", "$(BASE_URL)/home.htm?usr=$(user)&pwd=$(password)&persist=y")

    split(String(r.body), "$(COOKIE_NAME)=")[2]
end
