module WebSocketClient

#= The original DandelionWebSockets.jl package is licensed under the MIT "Expat" License: =#
#=  =#
#= > Copyright (c) 2016: Erik Edin. =#
#= > =#
#= > Permission is hereby granted, free of charge, to any person obtaining =#
#= > a copy of this software and associated documentation files (the =#
#= > "Software"), to deal in the Software without restriction, including =#
#= > without limitation the rights to use, copy, modify, merge, publish, =#
#= > distribute, sublicense, and/or sell copies of the Software, and to =#
#= > permit persons to whom the Software is furnished to do so, subject to =#
#= > the following conditions: =#
#= > =#
#= > The above copyright notice and this permission notice shall be =#
#= > included in all copies or substantial portions of the Software. =#
#= > =#
#= > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, =#
#= > EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF =#
#= > MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. =#
#= > IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY =#
#= > CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, =#
#= > TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE =#
#= > SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. =#

export AbstractWSClient,
       WSClient,
       stop,
       send_text,
       send_binary

export WebSocketHandler,
       on_text,
       on_binary,
       state_closed,
       state_closing,
       state_connecting,
       state_open,
       wsconnect

abstract type AbstractWSClient end

# This defines the public interface that the user should implement. These are callbacks called when
# events arrive from this WebSocket library.
abstract type WebSocketHandler end

"Handle a text frame."
on_text(t::WebSocketHandler, ::String) = nothing

"Handle a binary frame."
on_binary(t::WebSocketHandler, ::Vector{UInt8}) = nothing

"The WebSocket was closed."
state_closed(t::WebSocketHandler) = nothing

"The WebSocket is about to close."
state_closing(t::WebSocketHandler) = nothing

"The WebSocket is trying to connect."
state_connecting(t::WebSocketHandler) = nothing

"The WebSocket is open and ready to send and receive messages."
state_open(t::WebSocketHandler) = nothing

# include("core.jl")

# Core defines the core WebSocket types, such as a frame and opcodes.

# Description of a WebSocket frame from https://tools.ietf.org/html/rfc6455, chapter 5.2.
#
#      0                   1                   2                   3
#      0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1
#     +-+-+-+-+-------+-+-------------+-------------------------------+
#     |F|R|R|R| opcode|M| Payload len |    Extended payload length    |
#     |I|S|S|S|  (4)  |A|     (7)     |             (16/64)           |
#     |N|V|V|V|       |S|             |   (if payload len==126/127)   |
#     | |1|2|3|       |K|             |                               |
#     +-+-+-+-+-------+-+-------------+ - - - - - - - - - - - - - - - +
#     |     Extended payload length continued, if payload len == 127  |
#     + - - - - - - - - - - - - - - - +-------------------------------+
#     |                               |Masking-key, if MASK set to 1  |
#     +-------------------------------+-------------------------------+
#     | Masking-key (continued)       |          Payload Data         |
#     +-------------------------------- - - - - - - - - - - - - - - - +
#     :                     Payload Data continued ...                :
#     + - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +
#     |                     Payload Data continued ...                |
#     +---------------------------------------------------------------+

# TODO: Rethink what we're exporting here. Should we export anything? None of this
#       should be part of the interface to the user.
export Frame,
       Opcode,
       OPCODE_CONTINUATION,
       OPCODE_TEXT,
       OPCODE_BINARY,
       OPCODE_CLOSE,
       OPCODE_PING,
       OPCODE_PONG

import Base.==

# TODO: Documentation.

immutable Opcode
    op::UInt8
end

const OPCODE_CONTINUATION = Opcode(0)
const OPCODE_TEXT         = Opcode(1)
const OPCODE_BINARY       = Opcode(2)
const OPCODE_CLOSE        = Opcode(8)
const OPCODE_PING         = Opcode(9)
const OPCODE_PONG         = Opcode(10)

==(a::Opcode, b::Opcode) = a.op == b.op

immutable Frame
    fin::Bool
    rsv1::Bool
    rsv2::Bool
    rsv3::Bool
    opcode::Opcode
    ismasked::Bool
    len::UInt8 # Is actually 7 bits.
    extended_len::UInt64
    mask::Array{UInt8}
    payload::Array{UInt8}
end

Frame(fin::Bool, opcode::Opcode, ismasked::Bool, len::Int,
      extended_len::Int, mask::Vector{UInt8}, payload::Vector{UInt8}) =
  Frame(fin, false, false, false, opcode, ismasked, len, extended_len, mask, payload)

==(a::Frame, b::Frame) = a.fin == b.fin && a.rsv1 == b.rsv1 && a.rsv2 == b.rsv2 &&
  a.opcode == b.opcode && a.ismasked == b.ismasked && a.len == b.len && a.extended_len == b.extended_len &&
  a.mask == b.mask && a.payload == b.payload

function Base.read(s::IO, ::Type{Frame})
  x    = read(s, UInt8)
  fin  = x & 0b1000_0000 != 0
  rsv1 = x & 0b0100_0000 != 0
  rsv2 = x & 0b0010_0000 != 0
  rsv3 = x & 0b0001_0000 != 0
  op   = x & 0b0000_1111

  y        = read(s, UInt8)
  ismasked = y & 0b1000_0000 != 0
  len      = y & 0b0111_1111

  extended_len::UInt64 = 0
  if len == 126
    extended_len = ntoh(read(s, UInt16))
  elseif len == 127
    extended_len = ntoh(read(s, UInt64))
  end

  mask = Array{UInt8,1}()
  if ismasked
    mask = read(s, UInt8, 4)
  end

  payload_length = extended_len != 0 ? extended_len : len
  payload = read(s, UInt8, payload_length)

  Frame(fin, rsv1, rsv2, rsv3, Opcode(op), ismasked, len, extended_len, mask, payload)
end

function Base.write(s::IO, frame::Frame)
  x1 =
    UInt8(frame.fin)  << 7 |
    UInt8(frame.rsv1) << 6 |
    UInt8(frame.rsv2) << 5 |
    UInt8(frame.rsv3) << 4 |
    frame.opcode.op & 0b0000_1111

  x2 = UInt8(frame.ismasked) << 7 |
       frame.len & 0b0111_1111

  write(s, x1)
  write(s, x2)

  if frame.len == 126
    write(s, hton(UInt16(frame.extended_len)))
  elseif frame.len == 127
    write(s, hton(frame.extended_len))
  end

  if frame.ismasked
    write(s, frame.mask)
  end

  write(s, frame.payload)
end


# include("taskproxy.jl")

# A `TaskProxy` is a proxy object for another object. A set of predefined functions are called on
# the task proxy, and those functions and their arguments are sent via a channel to another
# coroutine, that performs the function calls on the target object.

abstract type TaskProxy end

const ProxyCall = Tuple{Function, Vector{Any}}

macro taskproxy(proxy_type::Symbol, abstract_type::Symbol, target_type::Symbol, functions...)

    proxy_functions = []
    # For each function in the macro arguments, make a function that collects the arguments and
    # sends the symbol and args on a channel.
    for fname in functions
        fexpr = :($fname(p::$proxy_type, args...) = put!(p.chan, ($fname, collect(args))))
        push!(proxy_functions, fexpr)
    end

    esc(
        quote
            # Define the proxy type, which contains the target object this acts as a proxy for, and
            # the channel that functions and arguments are sent over.
            # The target object can be unset at the beginning, and set with a call to `attach`
            # later on.
            type $proxy_type <: $abstract_type
                target::Nullable{$target_type}
                chan::Channel{ProxyCall}

                $(proxy_type)() = new(Nullable{$target_type}(), Channel{ProxyCall}(32))
                $(proxy_type)(target::$target_type) =
                    new(Nullable{$target_type}(target), Channel{ProxyCall}(32))
            end

            $(proxy_functions...)

            function do_proxy(p::$proxy_type, target::$target_type)
                for (f, args) in p.chan
                  try
                    f(target, args...)
                  catch ex
                    # task can fail when socket closed
                    if isa(ex, ArgumentError)
                      println("WARN: " * string(ex.msg))
                    else
                      println("WARN: " * string(ex))
                    end
                  end
                end
            end

            function start(p::$proxy_type)
                if isnull(p.target)
                    error("Target not set in proxy $(p). Call `attach` or set in constructor")
                end
                target = get(p.target)
                @async do_proxy(p, target)
            end

            function stop(h::$proxy_type)
                close(h.chan)
                h.target = Nullable{$target_type}()
            end

            is_set(p::$proxy_type) = !isnull(p.target)

            "Set the target object for an empty task proxy."
            function attach(p::$proxy_type, target::$target_type)
                !isnull(p.target) && error("Target already set")
                p.chan = Channel{ProxyCall}(32)
                p.target = Nullable{$target_type}(target)
            end
        end
    )
end


# include("glue_interface.jl")

# This file merely define abstract types used by the proxies that glue the different parts together.
export AbstractClientLogic

abstract type AbstractHandlerTaskProxy <: TaskProxy end
abstract type AbstractClientTaskProxy <: TaskProxy end
abstract type AbstractWriterTaskProxy <: TaskProxy end

abstract type AbstractClientLogic end
abstract type AbstractPinger end
abstract type AbstractPonger end


# include("network.jl")

import Base: read, write

# TODO: Documentation

"An exception thrown into a task in order to stop it."
type StopTaskException <: Exception end

abstract type AbstractServerReader end

"Reading from a network socket and placing the resulting frame on a channel."
immutable ServerReader <: AbstractServerReader
    s::IO
    task::Task
end

"Read frames from the network, until an exception is thrown in this task."
function do_reader(s::IO, logic::AbstractClientTaskProxy)
    try
        while true
            frame = read(s, Frame)
            # This is a proxy, so the actual `handle` call made on the logic object is done in a
            # separate coroutine.
            handle(logic, FrameFromServer(frame))
        end
    catch ex
        # TODO: Handle errors better.
    end
    handle(logic, SocketClosed())
end

function start_reader(s::IO, logic::AbstractClientTaskProxy)
    t = @schedule do_reader(s, logic)
    ServerReader(s, t)
end


function stop(t::ServerReader)
    try
        Base.throwto(t.task, StopTaskException())
    end
end

#
#
#

"""
TLSBufferedIO adapts a TLS socket so we can do byte I/O.

The stream returned by MbedTLS when using a TLS socket does not support the byte I/O used when
reading a frame. It only supports reading a chunk of data. This is a fake stream that buffers some
data and lets us do byte I/O.

Note: This should have been done by the BufferedStreams.jl package. However, I couldn't get it to
work with the MbedTLS stream, for reasons unknown. If we can investigate and fix that problem, then
we should really replace this type with a BufferedInputStream.
"""
immutable TLSBufferedIO <: IO
    tls_stream::IO
    buf::IOBuffer

    TLSBufferedIO(tls_stream::IO) = new(tls_stream, IOBuffer())
end

"Read all available data, and block until we have enough to fulfíll the next read."
function fill_buffer(s::TLSBufferedIO, n::Int)
    begin_ptr = mark(s.buf)
    while s.buf.size - begin_ptr < n
        write(s.buf, readavailable(s.tls_stream))
    end
    reset(s.buf)
end

function read(s::TLSBufferedIO, t::Type{UInt8})
    fill_buffer(s, sizeof(t))
    read(s.buf, t)
end

function read(s::TLSBufferedIO, t::Type{UInt16})
    fill_buffer(s, sizeof(t))
    read(s.buf, t)
end

function read(s::TLSBufferedIO, t::Type{UInt64})
    fill_buffer(s, sizeof(t))
    read(s.buf, t)
end

function read(s::TLSBufferedIO, t::Type{UInt8}, n::Int)
    fill_buffer(s, sizeof(t) * n)
    read(s.buf, t, n)
end

write(s::TLSBufferedIO, t::UInt8) = write(s.tls_stream, t)
write(s::TLSBufferedIO, t::UInt16) = write(s.tls_stream, t)
write(s::TLSBufferedIO, t::UInt64) = write(s.tls_stream, t)


# include("client_logic.jl")

# Client logic deals with handling control frames, user requesting to send frames, state.
# The ClientLogic type is defined below entirely synchronously. It takes input via the `handle()`
# function, which is defined for the different input types below. It performs internal logic and
# produces a call to its outbound interface.
#
# The outbound interface is composed of two abstract types `AbstractHandlerTaskProxy` and
# `AbstractWriterTaskProxy`. The concrete implementations will be `TaskProxy` objects, which will
# store the calls (function and arguments) and call it on a target object in another coroutine. This
# means that as long as the channels don't block, the logic will be performed concurrently with the
# the callbacks and writing to the network. This is important because the logic might have to respond
# to ping requests in a timely manner, which it might not be able to do if the callbacks block.
#
# For testing purposes the abstract outbound interface can be replaced with mock objects. This lets us
# test the logic of the WebSocket synchronously, without any asynchronicity or concurrency
# complicating things.

export ClientLogic

#
# These types define the input interface for the client logic.
#

"Abstract type for all commands sent to `ClientLogic`.

These commands are sent as arguments to the different `handle` functions on `ClientLogic`. Each
command represents an action on a WebSocket, such as sending a text frame, ping request, or closing
the connection."
abstract type ClientLogicInput end

"Send a text frame, sent to `ClientLogic`."
immutable SendTextFrame <: ClientLogicInput
	data::String
	# True if this is the final frame in the text message.
	isfinal::Bool
	# What WebSocket opcode should be used.
	opcode::Opcode
end

"Send a binary frame, sent to `ClientLogic`."
immutable SendBinaryFrame <: ClientLogicInput
	data::Array{UInt8, 1}
	# True if this is the final frame in the text message.
	isfinal::Bool
	# What WebSocket opcode should be used.
	opcode::Opcode
end

"Send a ping request to the server."
immutable ClientPingRequest  <: ClientLogicInput end

"A frame was received from the server."
immutable FrameFromServer <: ClientLogicInput
	frame::Frame
end

"A request to close the WebSocket."
immutable CloseRequest <: ClientLogicInput end

"Used when the underlying network socket was closed."
immutable SocketClosed <: ClientLogicInput end

"A pong reply was expected, but never received."
immutable PongMissed <: ClientLogicInput end

#
# ClientLogic
#

"Enum value for the different states a WebSocket can be in."
immutable SocketState
	v::Symbol
end

# We never send a `state_connecting` callback here, because that should be done when we make the
# HTTP upgrade.
const STATE_CONNECTING     = SocketState(:connecting)
# We send a `state_open` callback when the ClientLogic is created, when making the connection.
const STATE_OPEN           = SocketState(:open)
const STATE_CLOSING        = SocketState(:closing)
const STATE_CLOSING_SOCKET = SocketState(:closing_socket)
const STATE_CLOSED         = SocketState(:closed)

"Type for the logic of a client WebSocket."
type ClientLogic <: AbstractClientLogic
	# A WebSocket can be in a number of states. See the `STATE_*` constants.
	state::SocketState
	# The object to which callbacks should be made. This proxy will make the callbacks
	# asynchronously.
	handler::AbstractHandlerTaskProxy
	# A proxy for the stream where we write our frames.
	writer::AbstractWriterTaskProxy
	# Random number generation, used for masking frames.
	rng::AbstractRNG
	# Keeps track of when a pong is expected to be received from the server.
	ponger::AbstractPonger
	# Here we keep data collected when we get a message made up of multiple frames.
	buffer::Vector{UInt8}
	# This stores the type of the multiple frame message. This is the opcode of the first frame,
	# as the following frames have the OPCODE_CONTINUATION opcode.
	buffered_type::Opcode
	# This function cleans up the client when the connection is closed.
	client_cleanup::Function
end

ClientLogic(state::SocketState,
			handler::AbstractHandlerTaskProxy,
			writer::AbstractWriterTaskProxy,
                        rng::AbstractRNG,
                        ponger::AbstractPonger,
                        client_cleanup::Function) =
    ClientLogic(state, handler, writer, rng, ponger, Vector{UInt8}(), OPCODE_TEXT, client_cleanup)

"Send a frame to the other endpoint, using the supplied payload and opcode."
function send(logic::ClientLogic, isfinal::Bool, opcode::Opcode, payload::Vector{UInt8})
	# We can't send any frames in CLOSING or CLOSED.
	if logic.state != STATE_OPEN
		return
	end

	# Each frame is masked with four random bytes.
	mask    = rand(logic.rng, UInt8, 4)
	masking!(payload, mask)

	len::UInt64  = length(payload)
	extended_len = 0

	if 128 <= len <= 65536
		extended_len = len
		len = 126
	elseif 65536 + 1 <= len
		extended_len = len
		len = 127
	end

	# Create a Frame and write it to the underlying socket, via the `writer` proxy.
	frame = Frame(
		isfinal, false, false, false, opcode, true, len, extended_len, mask, payload)

	write(logic.writer, frame)
end

"Send a single text frame."
function handle(logic::ClientLogic, req::SendTextFrame)
	payload = Vector{UInt8}(req.data)
	send(logic, req.isfinal, req.opcode, payload)
end

"Send a single binary frame."
handle(logic::ClientLogic, req::SendBinaryFrame)   = send(logic, req.isfinal, req.opcode, req.data)

function handle(logic::ClientLogic, req::ClientPingRequest)
	if logic.state == STATE_OPEN
		ping_sent(logic.ponger)
		send(logic, true, OPCODE_PING, b"")
	end
end

function handle(logic::ClientLogic, ::PongMissed)
	logic.state = STATE_CLOSED
	state_closed(logic.handler)
	logic.client_cleanup()
end

"Handle a user request to close the WebSocket."
function handle(logic::ClientLogic, req::CloseRequest)
	logic.state = STATE_CLOSING

	# Send a close frame to the server
	mask = rand(logic.rng, UInt8, 4)
	frame = Frame(true, false, false, false, OPCODE_CLOSE, true, 0, 0,
		mask, b"")
	write(logic.writer, frame)

	state_closing(logic.handler)
end

"The underlying socket was closed. This is sent by the reader."
function handle(logic::ClientLogic, ::SocketClosed)
	logic.state = STATE_CLOSED
	state_closed(logic.handler)
	logic.client_cleanup()
end

"Handle a frame from the server."
function handle(logic::ClientLogic, req::FrameFromServer)
	if req.frame.opcode == OPCODE_CLOSE
		handle_close(logic, req.frame)
	elseif req.frame.opcode == OPCODE_PING
		handle_ping(logic, req.frame.payload)
	elseif req.frame.opcode == OPCODE_PONG
		handle_pong(logic, req.frame.payload)
	elseif req.frame.opcode == OPCODE_TEXT
		handle_text(logic, req.frame)
	elseif req.frame.opcode == OPCODE_BINARY
		handle_binary(logic, req.frame)
	elseif req.frame.opcode == OPCODE_CONTINUATION
		handle_continuation(logic, req.frame)
	end
end

#
# Internal handle functions
#

function handle_close(logic::ClientLogic, frame::Frame)
	# If the server initiates a closing handshake when we're in open, we should reply with a close
	# frame. If the client initiated the closing handshake then we'll be in STATE_CLOSING when the
	# reply comes, and we shouldn't send another close frame.
	send_close_reply = logic.state == STATE_OPEN
	logic.state = STATE_CLOSING_SOCKET
	if send_close_reply
		mask = rand(logic.rng, UInt8, 4)
		frame = Frame(true, false, false, false, OPCODE_CLOSE, true, frame.len, frame.extended_len,
			mask, frame.payload)
		write(logic.writer, frame)
		state_closing(logic.handler)
	end
end

function handle_ping(logic::ClientLogic, payload::Vector{UInt8})
	send(logic, true, OPCODE_PONG, payload)
end

function handle_pong(logic::ClientLogic, ::Vector{UInt8})
	pong_received(logic.ponger)
end

function handle_text(logic::ClientLogic, frame::Frame)
	if frame.fin
		on_text(logic.handler, String(frame.payload))
	else
		start_buffer(logic, frame.payload, OPCODE_TEXT)
	end
end

function handle_binary(logic::ClientLogic, frame::Frame)
	if frame.fin
		on_binary(logic.handler, frame.payload)
	else
		start_buffer(logic, frame.payload, OPCODE_BINARY)
	end
end

# TODO: What if we get a binary/text frame before we get a final continuation frame?
function handle_continuation(logic::ClientLogic, frame::Frame)
	buffer(logic, frame.payload)
	if frame.fin
		if logic.buffered_type == OPCODE_TEXT
			on_text(logic.handler, String(logic.buffer))
		elseif logic.buffered_type == OPCODE_BINARY
			on_binary(logic.handler, logic.buffer)
			logic.buffer = Vector{UInt8}()
		end
	end
end

function start_buffer(logic::ClientLogic, payload::Vector{UInt8}, opcode::Opcode)
	logic.buffered_type = opcode
	logic.buffer = copy(payload)
end

function buffer(logic::ClientLogic, payload::Vector{UInt8})
	append!(logic.buffer, payload)
end

#
# Utilities
#

function masking!(input::Vector{UInt8}, mask::Vector{UInt8})
	m = 1
	for i in 1:length(input)
		input[i] = input[i] ⊻ mask[(m - 1) % 4 + 1]
		m += 1
	end
end


# include("ping.jl")

export Pinger, stop,
       Ponger, pong_received, attach, ping_sent

type Pinger <: AbstractPinger
    timer::Nullable{Timer}
    interval::Float64

    function Pinger(interval::Float64)
        new(Nullable{Timer}(), interval)
    end
end

function attach(pinger::Pinger, logic::AbstractClientTaskProxy)
    send_ping = x -> handle(logic, ClientPingRequest())
    pinger.timer = Nullable{Timer}(Timer(send_ping, pinger.interval, pinger.interval))
end

function stop(p::Pinger)
    if !isnull(p.timer)
        close(get(p.timer))
        p.timer = Nullable{Timer}()
    end
end

type Ponger <: AbstractPonger
    timeout::Float64
    pong_missed::Function
    pongs_received::UInt64
    misses::Int
    current_misses::Int

    Ponger(timeout::Float64; misses::Int=1) = new(timeout, x -> nothing, 0, misses, 0)
end

function start_timer_(p::Ponger)
    p.timer = Nullable{Timer}(Timer(p.pong_missed, p.timeout, p.timeout))
end

attach(ponger::Ponger, logic::AbstractClientTaskProxy) =
    ponger.pong_missed = () -> handle(logic, PongMissed())

function pong_received(ponger::Ponger)
    ponger.pongs_received += 1
    ponger.current_misses = 0
end

function ping_sent(ponger::Ponger)
    pongs_received_at_send = ponger.pongs_received
    fun = x -> begin
        ponger.current_misses += 1
        if ponger.pongs_received == pongs_received_at_send && ponger.current_misses >= ponger.misses
            ponger.pong_missed()
            ponger.current_misses = 0
        end
    end
    Timer(fun, ponger.timeout)
end


# include("handshake.jl")

import Nettle
import Requests

"Keeps the result of a HTTP Upgrade attempt, when converting a HTTP connection to a WebSocket."
immutable HandshakeResult
    # The expected `Sec-WebSocket-Accept` value. See `validate`.
    expected_accept::String

    # The network stream opened by Requests, that we'll use for the WebSocket protocol.
    stream::IO

    # Response headers, keeping the WebSocket accept value, among others.
    headers::Dict{String,String}

    # When doing the HTTP upgrade, we might have read a part of the first WebSocket frame. This
    # contains that data.
    body::Vector{UInt8}
end

# Currently only used to dispatch a failed HTTP upgrade to another function.
immutable HandshakeFailure

end

"""
The WebSocket server is expected to reply with a computed value to prove that it's actually a
WebSocket server, and not another server duped into accepting this HTTP upgrade. This function
validates that the expected computed value is found in the response headers.
"""
function validate(handshake::HandshakeResult)
    normal_keys = collect(keys(handshake.headers))
    lower_keys = map(lowercase, normal_keys)
    accept_name_index = findfirst(lower_keys, "sec-websocket-accept")
    if accept_name_index == 0
        println("No Sec-WebSocket-Accept in $(handshake.headers)")
        return false
    end

    accept_name = normal_keys[accept_name_index]
    accept_value = handshake.headers[accept_name]

    is_valid = accept_value == handshake.expected_accept
    if !is_valid
        println("Expected accept value $(handshake.expected_accept) does not match actual $accept_value")
    end

    is_valid
end

"Create a random key that the server will use to compute its response."
function make_websocket_key(rng::AbstractRNG)
    ascii(base64encode(rand(rng, UInt8, 16)))
end

"Calculate the accept value, given the random key supplied by the client."
function calculate_accept(key::String)
    magic = "258EAFA5-E914-47DA-95CA-C5AB0DC85B11"
    h = Nettle.digest("sha1", key * magic)
    base64encode(h)
end

"Create headers used to upgrade the HTTP connection to a WebSocket connection."
function make_headers(key::String)
    headers = Dict(
        "Upgrade" => "websocket",
        "Connection" => "Upgrade",
        "Sec-WebSocket-Key" => key,
        "Sec-WebSocket-Version" => "13")
end

"Make a HTTP connection and upgrade it to a WebSocket connection."
function do_handshake(rng::AbstractRNG, uri::Requests.URI; do_request=Requests.do_stream_request)
    key = make_websocket_key(rng)
    expected_accept = calculate_accept(key)
    headers = make_headers(key)
    result = do_request(uri, ascii("GET"); headers=headers)

    stream = result.socket
    if uri.scheme == "https"
        stream = TLSBufferedIO(stream)
    end

    # TODO: Any body unintentionally read during the HTTP parsing is not returned, which means that
    #       if any such bytes were read, then we will not be able to correctly read the first frame.
    HandshakeResult(expected_accept, stream, result.response.headers, b"")
end

"Convert `ws://` or `wss://` URIs to 'http://` or `https://`."
function convert_ws_uri(uri::Requests.URI)
    u = replace(string(uri), r"^ws", "http")
    Requests.URI(u)
end


# include("client.jl")

import Requests: URI
import Base: show
using BufferedStreams

# These proxies glue the different coroutines together. For isntance, `ClientLogic` calls callback
# function such as `on_text` and `state_closing` on the proxy, which is then called on the callback
# object by another coroutine. This lets the logic run independently of the callbacks.
@taskproxy(HandlerTaskProxy, AbstractHandlerTaskProxy, WebSocketHandler,
    on_text, on_binary,
    state_connecting, state_open, state_closing, state_closed)

@taskproxy ClientLogicTaskProxy AbstractClientTaskProxy AbstractClientLogic handle
@taskproxy WriterTaskProxy AbstractWriterTaskProxy IO write

"""
A WebSocket client, used to connect to a server, and send messages.

Note: The keyword arguments in the constructor are primarily for testing.
"""
type WSClient <: AbstractWSClient
    # `writer` writes frames to the socket.
    writer::AbstractWriterTaskProxy
    # `handler_proxy` does the callbacks, in its own coroutine.
    handler_proxy::AbstractHandlerTaskProxy
    # `logic_proxy` forwards commands to the `ClientLogic` object, in its own coroutine.
    logic_proxy::AbstractClientTaskProxy
    # `reader` reads frames from the server.
    reader::Nullable{ServerReader}
    # `do_handshake` is a function that performs a HTTP Upgrade to a WebSocket connection.
    do_handshake::Function
    # `rng` is used for random generation of masks when sending frames.
    rng::AbstractRNG
    # `ponger` keeps track of when a pong response is expected from the server.
    ponger::AbstractPonger
    # `pinger` requests that the logic send ping frames to the server at regular intervals.
    pinger::AbstractPinger

    function WSClient(;
                      do_handshake=WebSocketClient.do_handshake,
                      rng::AbstractRNG=MersenneTwister(0),
                      writer::AbstractWriterTaskProxy=WriterTaskProxy(),
                      handler_proxy::AbstractHandlerTaskProxy=HandlerTaskProxy(),
                      logic_proxy::AbstractClientTaskProxy=ClientLogicTaskProxy(),
                      ponger::AbstractPonger=Ponger(3.0; misses=3),
                      pinger::AbstractPinger=Pinger(5.0))
        new(writer, handler_proxy, logic_proxy, Nullable{ServerReader}(), do_handshake, rng, ponger, pinger)
    end
end
show(io::IO, c::WSClient) =
    show(io, "WSClient($(c.handler_proxy), $(c.logic_proxy))")

"Validates a HTTP Upgrade response, and starts all tasks.

Note: As of right now the handshake is not validated, because the response headers aren't set here.
"
function connection_result_(client::WSClient, result::HandshakeResult, handler::WebSocketHandler)
    # Validation of a HTTP Upgrade to a WebSocket is done by checking the response headers for a key
    # which should contain a computed value.
    if !validate(result)
        println("Could not validate HTTP Upgrade")
        state_closed(handler)
        return false
    end

    # Each `TaskProxy` used here acts as a proxy for another object. When you call some predefined
    # functions on a proxy, it takes the function and arguments and puts them on a channel. A
    # coroutine takes there function/arguments from the channel and calls the same function, but on
    # the target object it acts as a proxy for. This is because we want some parts to work
    # concurrently with others
    # Calling `attach` on a proxy sets the target object.
    # Calling `start` on a proxy starts the coroutine that calls functions on that target.

    # For `writer` the target object is the IO stream for the WebSocket connection.
    attach(client.writer, result.stream)
    start(client.writer)

    # For `handler_proxy` the target is the handler object on which callbacks should be made.
    attach(client.handler_proxy, handler)
    start(client.handler_proxy)

    # Note: This doesn't directly call the `state_open` callback on the handler, but rather enqueues
    # the function call, so that the `handler_proxy` coroutine will make the actual callback.
    state_open(client.handler_proxy)

    # This function stops all the task proxies, effectively cleaning up the WSClient. This is
    # necessary when one wants to reconnect.
    cleanup = () -> begin
        stop(client.writer)
        stop(client.handler_proxy)
        stop(client.logic_proxy)
        stop(client.pinger)
        if !isnull(client.reader)
            stop(get(client.reader))
        end
    end

    # `ClientLogic` starts in the `STATE_OPEN` state, because it isn't responsible for making
    # connections. The target object for `logic_proxy` is the `ClientLogic` object created here.
    logic = ClientLogic(STATE_OPEN, client.handler_proxy, client.writer, client.rng, client.ponger,
                        cleanup)
    attach(client.logic_proxy, logic)
    start(client.logic_proxy)

    # `Ponger` requires a logic object it can alert when a pong request hasn't been received within
    # the expected time frame. This attaches that logic object to the ponger.
    attach(client.ponger, client.logic_proxy)

    # `Pinger` sends ping requests at regular intervals.
    attach(client.pinger, client.logic_proxy)

    # The target for `reader` is the same stream we're writing to.
    client.reader = Nullable{ServerReader}(
        start_reader(result.stream, client.logic_proxy))
    true
end

"The HTTP Upgrade failed, for whatever reason."
function connection_result_(client::WSClient, result::HandshakeFailure, handler::WebSocketHandler)
    # Calling `state_closed` here, because that's where we're expected to attempt to reconnect.
    # Note: We call `state_closed` directly, not using the proxy, because the proxy hasn't been
    # started yet.
    state_closed(handler)
    false
end

"Connect the client to a WebSocket server at `uri`, and use `handler` for the callbacks."
function wsconnect(client::WSClient, uri::URI, handler::WebSocketHandler)
    # Note: Calling the `state_connecting` callback directly, because the `handler_proxy` hasn't
    # been started yet.
    state_connecting(handler)

    # This converts from `ws://` or `wss://` to `http://` or `https://`, because that's what
    # Requests.jl expects.
    new_uri = convert_ws_uri(uri)

    # This makes a HTTP request to the URI and attempts to upgrade the connection to the WebSocket
    # protocol.
    handshake_result = client.do_handshake(client.rng, new_uri)

    connection_result_(client, handshake_result, handler)
end

"Close the WebSocket connection."
stop(c::WSClient) = handle(c.logic_proxy, CloseRequest())

"Send a single text frame."
# make a copy of the string before passing to SendTextFrame because the input string will be modified by the method
send_text(c::WSClient, s::String) =
  handle(c.logic_proxy, SendTextFrame(s * "", true, OPCODE_TEXT))

"Send a single binary frame."
# make a copy of the byte array before passing to SendTextFrame because the byte array will be modified by the method
send_binary(c::WSClient, data::Vector{UInt8}) =
    handle(c.logic_proxy, SendBinaryFrame(copy(data), true, OPCODE_BINARY))

end # module
