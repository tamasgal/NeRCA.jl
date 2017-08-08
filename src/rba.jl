import Requests: URI
using JSON

#const RBA_URI = URI("ws://maroff1.in2p3.fr:8088/message")
const RBA_URI = URI("ws://localhost:8088/message")


type RBAHandler <: WebSocketClient.WebSocketHandler
    client::WebSocketClient.WSClient
    stop_channel::Queue{Any}
end


function rba(token::AbstractString, hits::Vector{CalibratedHit})
    raw = Dict("pos" => [h.pos for h in hits],
               "time" => [h.t for h in hits],
               "tot" => [Int(h.tot) for h in hits])

    message = Dict("token" => token,
   		   "kind" => "event",
		   "data" => Dict("hits" => raw))

    handler = RBAHandler(WebSocketClient.WSClient(), Queue(Any))

    WebSocketClient.state_connecting(::RBAHandler) = println("Connecting to RBA server...")
    WebSocketClient.state_open(handler::RBAHandler) = println("connection to RBA successful.")
    WebSocketClient.state_closing(::RBAHandler) = println("Closing connection to RBA...")
    WebSocketClient.state_closed(handler::RBAHandler) = println("connection closed.")


    WebSocketClient.wsconnect(handler.client, RBA_URI, handler)
    WebSocketClient.send_text(handler.client, JSON.json(message))
    WebSocketClient.stop(handler.client)
end
