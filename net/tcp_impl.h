/*
 * TcpClient
 * @author tweber
 * @date aug-2014
 * @license GPLv2 or GPLv3
 * @desc based on chat client by C. M. Kohlhoff:
 * @desc	www.boost.org/doc/libs/1_56_0/doc/html/boost_asio/example/cpp11/chat/chat_client.cpp
 */

#ifndef __TL_TCP_IMPL_H__
#define __TL_TCP_IMPL_H__

#include "tcp.h"
#include "../log/log.h"
#include <boost/tokenizer.hpp>

namespace tl {

template<class t_ch, class t_str>
bool TcpTxtClient_gen<t_ch, t_str>::get_cmd_tokens(const t_str& str, const t_str& strDelim,
	std::vector<t_str>& vecStr, t_str& strRemainder)
{
	boost::char_separator<t_ch> delim(strDelim.c_str(), "", boost::keep_empty_tokens);
	boost::tokenizer<boost::char_separator<t_ch>> tok(str, delim);

	for(const t_str& strTok : tok)
	{
		//tl::log_debug("tok: ", strTok);
		vecStr.push_back(strTok);
	}

	if(vecStr.size()<=1)
		return false;

	// keep_empty_tokens leads to an empty last element if no remaining string is left
	if(*vecStr.rbegin() == "")
	{
		strRemainder = "";
		vecStr.pop_back();
	}
	else
	{
		strRemainder = *vecStr.rbegin();
		vecStr.pop_back();
	}

	//tl::log_debug("toks: ", vecStr.size());
	//tl::log_debug("remainder: ", strRemainder);
	return true;
}


template<class t_ch, class t_str>
TcpTxtClient_gen<t_ch, t_str>::TcpTxtClient_gen() : m_listWriteBuffer(1024)
{}

template<class t_ch, class t_str>
TcpTxtClient_gen<t_ch, t_str>::~TcpTxtClient_gen()
{
	disconnect();

	m_sigRecv.disconnect_all_slots();
	m_sigDisconn.disconnect_all_slots();
	m_sigConn.disconnect_all_slots();

	const t_str* pstr = nullptr;
	while(m_listWriteBuffer.pop(pstr))
		if(pstr) delete pstr;
}

template<class t_ch, class t_str>
bool TcpTxtClient_gen<t_ch, t_str>::connect(const t_str& strHost, const t_str& strService)
{
	m_strHost = strHost;
	m_strService = strService;

	try
	{
		disconnect();

		m_pservice = new asio::io_service;
		m_psock = new ip::tcp::socket(*m_pservice);

		ip::tcp::resolver res(*m_pservice);
		ip::tcp::resolver::iterator iter = res.resolve({strHost, strService});

		asio::async_connect(*m_psock, iter,
		[this](const sys::error_code& err, ip::tcp::resolver::iterator)
		{
			if(!err)
				read_loop();

			m_sigConn(m_strHost, m_strService);
		});

		//tl::log_debug("Starting service thread.");
		m_pthread = new std::thread([this]()
		{
			try
			{
				m_pservice->run();
			}
			catch(const std::exception& ex)
			{
				log_err("TCP service thread exited with error: ", ex.what(), ".");
				m_pthread = nullptr;
				disconnect(1);
			}
		});
	}
	catch(const std::exception& ex)
	{
		log_err(ex.what());
		return 0;
	}

	return 1;
}

template<class t_ch, class t_str>
void TcpTxtClient_gen<t_ch, t_str>::disconnect(bool bAlwaysSendSignal)
{
	const bool bConnected = is_connected();
	if(bConnected)
	{
		m_pservice->stop();
		m_psock->close();
	}

	if(m_psock) { delete m_psock; m_psock = 0; }
	if(m_pthread)
	{
		m_pthread->join();
		delete m_pthread;
		m_pthread = 0;
	}
	if(m_pservice) { delete m_pservice; m_pservice = 0; }

	if(bConnected || bAlwaysSendSignal)
	{
		m_sigDisconn(m_strHost, m_strService);
		m_strHost = m_strService = "";
	}
}

template<class t_ch, class t_str>
bool TcpTxtClient_gen<t_ch, t_str>::is_connected()
{
	if(!m_psock) return 0;
	return m_psock->is_open();
}

template<class t_ch, class t_str>
void TcpTxtClient_gen<t_ch, t_str>::wait()
{
	if(m_pthread)
		m_pthread->join();
}

template<class t_ch, class t_str>
void TcpTxtClient_gen<t_ch, t_str>::write(const t_str& str)
{
	//tl::log_debug("write: ", str);
	m_listWriteBuffer.push(new t_str(str));
	m_pservice->post([&](){ flush_write(); });
}

template<class t_ch, class t_str>
void TcpTxtClient_gen<t_ch, t_str>::flush_write()
{
	const t_str* pstr = nullptr;
	if(m_listWriteBuffer.empty()) return;
	if(!m_listWriteBuffer.pop(pstr)) return;
	if(!pstr) return;
	const t_str& str = *pstr;

	//tl::log_debug("str = ", str);
	asio::async_write(*m_psock, asio::buffer(str.data(), str.length()),
	[this, pstr](const sys::error_code& err, std::size_t len)
	{
		if(pstr) delete pstr;

		if(err)
		{
			disconnect();
			return;
		}

		if(!m_listWriteBuffer.empty())
			flush_write();
	});
}

template<class t_ch, class t_str>
void TcpTxtClient_gen<t_ch, t_str>::read_loop()
{
	static const std::size_t iBufLen = 512;
	static t_ch pcBuf[iBufLen];

	//tl::log_debug("In read loop.");
	asio::async_read(*m_psock, asio::buffer(pcBuf, iBufLen), asio::transfer_at_least(1),
	[this](const sys::error_code& err, std::size_t len)
	{
		if(err)
		{
			log_err("TCP read error.");
			disconnect();
			return;
		}

		t_str strCurMsg(pcBuf, len);
		m_strReadBuffer.append(strCurMsg);

		//tl::log_debug("read buffer: ", m_strReadBuffer);
		std::vector<t_str> vecCmds;
		if(get_cmd_tokens(m_strReadBuffer, m_strCmdDelim, vecCmds, m_strReadBuffer))
		{
			//tl::log_debug("remainder: ", m_strReadBuffer);
			for(const t_str& strCmd : vecCmds)
			{
				//tl::log_debug("Read: ", strCmd);
				m_sigRecv(strCmd);
			}
		}

		read_loop();
	});
  }


// --------------------------------------------------------------------------------
// Signals
template<class t_ch, class t_str>
void TcpTxtClient_gen<t_ch, t_str>::add_receiver(const typename t_sigRecv::slot_type& conn)
{
	m_sigRecv.connect(conn);
}

template<class t_ch, class t_str>
void TcpTxtClient_gen<t_ch, t_str>::add_disconnect(const typename t_sigDisconn::slot_type& disconn)
{
	m_sigDisconn.connect(disconn);
}

template<class t_ch, class t_str>
void TcpTxtClient_gen<t_ch, t_str>::add_connect(const typename t_sigConn::slot_type& conn)
{
	m_sigConn.connect(conn);
}
// --------------------------------------------------------------------------------

}

#endif
