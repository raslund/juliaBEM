module statusPrints

    export message, startup


    function message(strArr; lMax = 60, caps="=",j="c")
        println(caps^lMax)
        for s in 1:length(strArr)
            if j == "c"
                padLength = Int(floor((lMax-length(strArr[s]))/2))
            elseif j == "l"
                padLength = 1
            end
            println(" "^padLength*strArr[s])
        end
        println(caps^lMax)
    end


    function startup()
        message(["JuliaBEM v0.1","RSL 26-03-2020"])
    end

end
