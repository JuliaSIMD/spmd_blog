using Dates

function hfun_bar(vname)
    val = Meta.parse(vname[1])
    return round(sqrt(val), digits=2)
end

function hfun_m1fill(vname)
    var = vname[1]
    return pagevar("index", var)
end

function lx_baz(com, _)
    # keep this first line
    brace_content = Franklin.content(com.braces[1]) # input string
    # do whatever you want here
    return uppercase(brace_content)
end

"""
    {{blogposts}}
Plug in the list of blog posts contained in the `/posts/` folder.
"""
function hfun_blogposts()
    io = IOBuffer()
    posts = readdir("posts")
    post_data = map(posts) do post
        ps  = splitext(post)[1]
        url = "/posts/$ps/"
        surl = string(strip(url, '/'))::String
        title = string(pagevar(surl, :title))::String
        d = pagevar(surl, :date)
        if d isa Date
           return d, title, url
        else
           throw("Date not found for: $surl\nfound pubdate = $d")
        end
    end
    sort!(post_data, by=first, rev=true)

    prev_yr = nothing
    for post in post_data
        pubdate, title, url = post
        yr = Dates.year(pubdate)
        if yr != prev_yr
            write(io, "\\section{$yr}\n")
        end
        prev_yr = yr
        write(io, "- `$pubdate` -- [$title]($url)\n")
    end
    return Franklin.fd2html(String(take!(io)), internal=true)
end

function lx_section(com, _)
    content = Franklin.content(com.braces[1])
    return """~~~<h3 class="b normal">~~~$content~~~</h3>~~~"""
end
